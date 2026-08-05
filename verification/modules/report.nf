process REPORT {
    tag 'aggregate'
    publishDir params.outdir, mode: 'copy'

    input:
    path comparison_files, stageAs: 'comparison???.json'

    output:
    path 'verification.json', emit: verification
    path 'report.csv', emit: report

    script:
    """
    python3 - comparison???.json <<'PY'
    import csv
    import json
    import pathlib
    import sys

    comparisons = [json.loads(pathlib.Path(name).read_text()) for name in sys.argv[1:]]
    comparisons.sort(key=lambda item: (item['lane'], item['case_id']))
    verification = {
        'schema_version': 1,
        'ok': all(item['ok'] for item in comparisons),
        'comparison_count': len(comparisons),
        'comparisons': comparisons,
    }
    pathlib.Path('verification.json').write_text(json.dumps(verification, indent=2, sort_keys=True) + '\\n')
    with open('report.csv', 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['lane', 'case_id', 'ok', 'difference_count'])
        for item in comparisons:
            writer.writerow([item['lane'], item['case_id'], item['ok'], len(item['differences'])])
    PY
    """
}
