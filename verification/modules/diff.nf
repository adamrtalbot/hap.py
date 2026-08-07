process DIFF_OUTPUTS {
    tag { "${case_name}:${meta.id}" }
    publishDir { "${params.outdir}/${case_name}/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), val(case_name), val(prefix), path(legacy_files, stageAs: 'legacy/*'), path(rust_files, stageAs: 'rust/*')

    output:
    tuple val(meta), val(case_name), path('comparison.json'), emit: comparison
    path 'legacy.jsonl.gz', emit: legacy_stream
    path 'rust.jsonl.gz', emit: rust_stream

    script:
    """
    python3 - '${case_name}' '${meta.id}' '${prefix}' <<'PY'
    import csv
    import gzip
    import hashlib
    import io
    import json
    import pathlib
    import re
    import subprocess
    import sys

    lane, case_id, prefix = sys.argv[1:]
    observed = {}
    differences = []

    def object_without_duplicate_keys(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError('duplicate JSON key: %s' % key)
            result[key] = value
        return result

    def content(path):
        value = path.read_bytes()
        return gzip.decompress(value) if path.name.endswith('.gz') else value

    volatile_json_pointers = {
        '/version',
        '/runInfo',
        '/timestamp',
        '/dist',
        '/environment',
        '/mac_ver',
        '/python_implementation',
        '/python_prefix',
        '/python_version',
        '/uname',
        '/metadata/required/version',
        '/metadata/required/description',
        '/final_args/engine_vcfeval_template',
    }
    volatile_csv_columns = {'sompyversion', 'sompycmd'}
    volatile_vcf_header = re.compile(
        r'^(?:##bcftools_[^=]*(?:Command|Version)=|##fileDate=|##source=)'
    )

    def normalize_json(value, pointer=''):
        if isinstance(value, dict):
            normalized = {}
            for key, child in value.items():
                child_pointer = pointer + '/' + key.replace('~', '~0').replace('/', '~1')
                if child_pointer not in volatile_json_pointers:
                    normalized[key] = normalize_json(child, child_pointer)
            return normalized
        if isinstance(value, list):
            return [normalize_json(child, pointer + '/' + str(index)) for index, child in enumerate(value)]
        return value

    def json_difference(legacy, rust, pointer=''):
        if type(legacy) is not type(rust):
            return pointer, legacy, rust, 'JSON types differ'
        if isinstance(legacy, dict):
            if set(legacy) != set(rust):
                return pointer, sorted(legacy), sorted(rust), 'JSON keys differ'
            for key in sorted(legacy):
                difference = json_difference(legacy[key], rust[key], pointer + '/' + key.replace('~', '~0').replace('/', '~1'))
                if difference:
                    return difference
        elif isinstance(legacy, list):
            if len(legacy) != len(rust):
                return pointer, len(legacy), len(rust), 'JSON array lengths differ'
            for index, (legacy_item, rust_item) in enumerate(zip(legacy, rust)):
                difference = json_difference(legacy_item, rust_item, pointer + '/' + str(index))
                if difference:
                    return difference
        elif legacy != rust:
            return pointer, legacy, rust, 'JSON values differ'
        return None

    def text_difference(legacy, rust):
        legacy_lines = legacy.decode('utf-8').splitlines()
        rust_lines = rust.decode('utf-8').splitlines()
        for index, (legacy_line, rust_line) in enumerate(zip(legacy_lines, rust_lines), start=1):
            if legacy_line != rust_line:
                return '/line/%s' % index, legacy_line, rust_line
        if len(legacy_lines) != len(rust_lines):
            index = min(len(legacy_lines), len(rust_lines)) + 1
            return '/line/%s' % index, legacy_lines[index - 1] if index <= len(legacy_lines) else None, rust_lines[index - 1] if index <= len(rust_lines) else None
        return None

    def normalize_csv(value):
        rows = list(csv.reader(value.decode('utf-8').splitlines()))
        if not rows:
            return value
        indexes = [index for index, header in enumerate(rows[0]) if header not in volatile_csv_columns]
        normalized = []
        for row in rows:
            normalized.append([row[index] if index < len(row) else '' for index in indexes])
        output = []
        for row in normalized:
            stream = io.StringIO()
            csv.writer(stream, lineterminator='').writerow(row)
            output.append(stream.getvalue())
        return ('\\n'.join(output) + '\\n').encode('utf-8')

    def normalize_vcf(value):
        lines = value.decode('utf-8').splitlines()
        normalized = []
        for line in lines:
            if volatile_vcf_header.match(line) or line.startswith('##CL='):
                continue
            normalized.append(line)
        return ('\\n'.join(normalized) + '\\n').encode('utf-8') if normalized else b''

    def decode_bcf(path):
        result = subprocess.run(
            ['bcftools', 'view', '--no-version', '-Ov', str(path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return result.stdout

    def index_data_path(index_path):
        if index_path.name.endswith('.tbi'):
            return index_path.with_name(index_path.name[:-4])
        if index_path.name.endswith('.csi'):
            return index_path.with_name(index_path.name[:-4])
        raise ValueError('unsupported index extension: %s' % index_path.name)

    def validate_index(index_path):
        data_path = index_data_path(index_path)
        if not data_path.is_file():
            raise ValueError('index companion is missing: %s' % data_path.name)
        subprocess.run(
            ['bcftools', 'index', '--stats', str(data_path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def artifact_record(path):
        value = content(path)
        return {
            'name': path.name,
            'bytes': len(value),
            'sha256': hashlib.sha256(value).hexdigest(),
        }
    for side in ('legacy', 'rust'):
        names = sorted(path.name for path in pathlib.Path(side).glob(prefix + '*'))
        observed[side] = names
    if observed['legacy'] != observed['rust']:
        differences.append({
            'kind': 'artifact_set',
            'expected': observed['legacy'],
            'actual': observed['rust'],
        })
    for artifact in sorted(set(observed['legacy']).intersection(observed['rust'])):
        legacy_path = pathlib.Path('legacy') / artifact
        rust_path = pathlib.Path('rust') / artifact
        if not legacy_path.is_file() or not rust_path.is_file():
            continue
        if artifact.endswith(('.tbi', '.csi')):
            difference = {
                'artifact': artifact,
                'legacy_sha256': hashlib.sha256(legacy_path.read_bytes()).hexdigest(),
                'rust_sha256': hashlib.sha256(rust_path.read_bytes()).hexdigest(),
            }
            try:
                validate_index(legacy_path)
                validate_index(rust_path)
                continue
            except Exception as error:
                difference.update(kind='invalid', location='/', expected=None, actual=None, reason=str(error))
                differences.append(difference)
                continue
        legacy_content = content(legacy_path)
        rust_content = content(rust_path)
        if legacy_content == rust_content and not artifact.endswith('.bcf'):
            continue
        difference = {
            'artifact': artifact,
            'legacy_sha256': hashlib.sha256(legacy_content).hexdigest(),
            'rust_sha256': hashlib.sha256(rust_content).hexdigest(),
        }
        try:
            if artifact.endswith(('.json', '.json.gz')):
                legacy_json = json.loads(legacy_content.decode('utf-8'), object_pairs_hook=object_without_duplicate_keys)
                rust_json = json.loads(rust_content.decode('utf-8'), object_pairs_hook=object_without_duplicate_keys)
                json_diff = json_difference(normalize_json(legacy_json), normalize_json(rust_json))
                if json_diff is None:
                    continue
                pointer, expected_value, actual_value, reason = json_diff
                difference.update(kind='json', location=pointer, expected=expected_value, actual=actual_value, reason=reason)
            elif artifact.endswith('.bcf'):
                text_diff = text_difference(
                    normalize_vcf(decode_bcf(legacy_path)),
                    normalize_vcf(decode_bcf(rust_path)),
                )
                if text_diff is None:
                    continue
                location, expected_value, actual_value = text_diff
                difference.update(kind='bcf', location=location, expected=expected_value, actual=actual_value, reason='ordered BCF content differs')
            elif artifact.endswith(('.csv', '.csv.gz', '.txt', '.vcf', '.vcf.gz', '.bed', '.bed.gz', '.fai')):
                if artifact.endswith(('.csv', '.csv.gz')):
                    legacy_content = normalize_csv(legacy_content)
                    rust_content = normalize_csv(rust_content)
                elif artifact.endswith(('.vcf', '.vcf.gz')):
                    legacy_content = normalize_vcf(legacy_content)
                    rust_content = normalize_vcf(rust_content)
                text_diff = text_difference(legacy_content, rust_content)
                if text_diff is None:
                    continue
                location, expected_value, actual_value = text_diff
                difference.update(kind='text', location=location, expected=expected_value, actual=actual_value, reason='ordered text differs')
            else:
                difference.update(kind='binary', location='/sha256', expected=difference['legacy_sha256'], actual=difference['rust_sha256'], reason='binary content differs')
        except Exception as error:
            difference.update(kind='invalid', location='/', expected=None, actual=None, reason=str(error))
        differences.append(difference)
    for filename, side, directory in (('legacy.jsonl.gz', 'legacy', 'legacy'), ('rust.jsonl.gz', 'rust', 'rust')):
        payload = {
            'schema_version': 1,
            'side': side,
            'artifacts': [artifact_record(pathlib.Path(directory) / name) for name in observed[side]],
        }
        with gzip.open(filename, 'wt', encoding='utf-8') as handle:
            handle.write(json.dumps(payload, sort_keys=True, separators=(',', ':')) + '\\n')
    pathlib.Path('comparison.json').write_text(json.dumps({
        'schema_version': 1,
        'lane': lane,
        'case_id': case_id,
        'ok': not differences,
        'legacy_artifacts': observed['legacy'],
        'observed_artifacts': observed,
        'differences': differences,
    }, indent=2, sort_keys=True) + '\\n')
    PY
    """
}
