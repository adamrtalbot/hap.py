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
    from collections import Counter
    import gzip
    import hashlib
    import io
    import json
    import math
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
    roc_key_columns = ('Type', 'Subtype', 'Subset', 'Filter', 'Genotype', 'QQ.Field', 'QQ')
    roc_series_columns = roc_key_columns[:-1]
    roc_count_column = re.compile(
        r'(?:FP[.](?:gt|al)|(?:TRUTH|QUERY)[.](?:TOTAL|TP|FN|FP|UNK)(?:[.](?:ti|tv|het|homalt))?)'
    )

    def normalize_metrics_table(value):
        if value.get('type') != 'Table' or not isinstance(value.get('data'), list):
            return value
        columns = [column for column in value['data'] if isinstance(column, dict) and isinstance(column.get('values'), list)]
        if not columns:
            return value
        if str(value.get('id', '')).startswith('roc.'):
            for column in columns:
                column['values'] = []
            return value
        row_count = len(columns[0]['values'])
        synthetic_index = next((column for column in columns if column.get('id') == 'types'), None)
        if synthetic_index is not None and all(len(column['values']) == row_count for column in columns):
            synthetic_index['values'] = list(range(row_count))
        return value

    def normalize_json(value, pointer=''):
        if isinstance(value, dict):
            normalized = {}
            for key, child in value.items():
                child_pointer = pointer + '/' + key.replace('~', '~0').replace('/', '~1')
                if child_pointer not in volatile_json_pointers:
                    normalized[key] = normalize_json(child, child_pointer)
            return normalize_metrics_table(normalized)
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

    def normalized_csv_table(value):
        rows = list(csv.reader(value.decode('utf-8').splitlines()))
        if not rows:
            return (), []
        indexes = [index for index, header in enumerate(rows[0]) if header not in volatile_csv_columns]
        normalized = [
            tuple(row[index] if index < len(row) else '' for index in indexes)
            for row in rows
        ]
        return normalized[0], normalized[1:]

    def normalize_csv(value):
        header, rows = normalized_csv_table(value)
        stream = io.StringIO()
        csv.writer(stream, lineterminator='\\n').writerows([header] + rows if header else [])
        return stream.getvalue().encode('utf-8')

    def count_record(row, counts):
        return {
            'row': list(row) if row is not None else None,
            'count': counts[row] if row is not None else 0,
        }

    def multiplicity_difference(expected_row, actual_row, expected_counts, actual_counts):
        return '/rows', count_record(expected_row, expected_counts), count_record(actual_row, actual_counts), 'CSV row multiplicities differ'

    def row_multiplicity_difference(expected_rows, actual_rows):
        expected_counts = Counter(expected_rows)
        actual_counts = Counter(actual_rows)
        row = next((
            row for row in sorted(set(expected_counts).union(actual_counts))
            if expected_counts[row] != actual_counts[row]
        ), None)
        return multiplicity_difference(row, row, expected_counts, actual_counts) if row is not None else None

    def finite_number(value):
        try:
            return math.isfinite(float(value))
        except (TypeError, ValueError):
            return False

    def nearest_thresholds(extra, legacy_rows, indexes):
        qq_index = indexes['QQ']
        series_indexes = [indexes[name] for name in roc_series_columns]
        series = tuple(extra[index] for index in series_indexes)
        candidates = [
            row for row in legacy_rows
            if tuple(row[index] for index in series_indexes) == series
        ]
        threshold = float(extra[qq_index])
        lower = max(
            (row for row in candidates if float(row[qq_index]) < threshold),
            key=lambda row: float(row[qq_index]),
            default=None,
        )
        upper = min(
            (row for row in candidates if float(row[qq_index]) > threshold),
            key=lambda row: float(row[qq_index]),
            default=None,
        )
        return lower, upper

    def additional_threshold_difference(extra, rust_counts, legacy_rows, indexes):
        actual = count_record(extra, rust_counts)
        type_index = indexes['Type']
        qq_index = indexes['QQ']
        if extra[type_index] not in {'SNP', 'INDEL'} or not finite_number(extra[qq_index]):
            return '/rows', None, actual, 'unexpected non-threshold CSV row'
        lower_row, upper_row = nearest_thresholds(extra, legacy_rows, indexes)
        if lower_row is None or upper_row is None:
            return '/rows', None, actual, 'additional threshold is not bounded by legacy ROC points'
        count_indexes = [index for name, index in indexes.items() if roc_count_column.fullmatch(name)]
        for index in count_indexes:
            neighbours = (lower_row[index], extra[index], upper_row[index])
            if not all(finite_number(value) for value in neighbours):
                continue
            lower_value, value, upper_value = map(float, neighbours)
            if not min(lower_value, upper_value) <= value <= max(lower_value, upper_value):
                return '/rows', {
                    'lower': list(lower_row),
                    'upper': list(upper_row),
                }, actual, 'additional threshold breaks ROC count monotonicity'
        return None

    def csv_multiset_difference(legacy, rust):
        legacy_header, legacy_rows = normalized_csv_table(legacy)
        rust_header, rust_rows = normalized_csv_table(rust)
        if legacy_header != rust_header:
            return '/header', list(legacy_header), list(rust_header), 'CSV headers differ'
        indexes = {name: index for index, name in enumerate(legacy_header)}
        if not set(roc_key_columns).issubset(indexes):
            return row_multiplicity_difference(legacy_rows, rust_rows)

        type_index = indexes['Type']
        qq_index = indexes['QQ']
        comparable_legacy = [row for row in legacy_rows if row[type_index]]
        legacy_counts = Counter(comparable_legacy)
        rust_counts = Counter(rust_rows)
        differing_row = next((
            row for row in sorted(legacy_counts)
            if legacy_counts[row] != rust_counts[row]
        ), None)
        if differing_row is not None:
            semantic_indexes = [indexes[name] for name in roc_key_columns]
            replacement = next((
                row for row in rust_rows
                if all(row[index] == differing_row[index] for index in semantic_indexes)
            ), None)
            return multiplicity_difference(differing_row, replacement, legacy_counts, rust_counts)

        additional_rows = sorted(set(rust_counts) - set(legacy_counts))
        legacy_thresholds = [
            row for row in comparable_legacy
            if row[type_index] in {'SNP', 'INDEL'} and finite_number(row[qq_index])
        ]
        for extra in additional_rows:
            if rust_counts[extra] != 1:
                return multiplicity_difference(extra, extra, legacy_counts, rust_counts)
            difference = additional_threshold_difference(extra, rust_counts, legacy_thresholds, indexes)
            if difference:
                return difference
        return None

    def is_roc_csv(artifact):
        return artifact.startswith(prefix + '.roc.') and artifact.endswith(('.csv', '.csv.gz'))

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
            elif is_roc_csv(artifact):
                csv_diff = csv_multiset_difference(legacy_content, rust_content)
                if csv_diff is None:
                    continue
                location, expected_value, actual_value, reason = csv_diff
                difference.update(kind='csv', location=location, expected=expected_value, actual=actual_value, reason=reason)
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
