#!/usr/bin/env bash
set -euo pipefail

readonly nextflow_version="26.04.6"
readonly nextflow_sha256="61a755edbed743cfbb568f3a6c67af68481a2f6a4d6dffcc4295e51318968281"
readonly nextflow_jar_sha256="2ca0251ae2d749317d9fbe5fe191a1616b5f44b608224268924c71b32f5ed9e2"
readonly nf_test_version="0.9.5"
readonly nf_test_sha256="b7679eb90cdc9642bfa89e9634db02ce6e699de53ad35ef0e2f4847634fc1641"
readonly -a curl_options=(
    --fail
    --location
    --proto '=https'
    --tlsv1.2
    --retry 4
    --retry-all-errors
    --retry-delay 2
    --connect-timeout 15
    --max-time 300
)

readonly destination="${1:?usage: install-verification-tools.sh DESTINATION}"
readonly nextflow_home="${NXF_HOME:?NXF_HOME must identify the isolated Nextflow cache}"
readonly temporary="$(mktemp -d)"
trap 'rm -rf "$temporary"' EXIT

mkdir -p "$destination" "$nextflow_home/framework/$nextflow_version"

curl "${curl_options[@]}" \
    "https://www.nextflow.io/releases/v${nextflow_version}/nextflow" \
    --output "$temporary/nextflow"
curl "${curl_options[@]}" \
    "https://www.nextflow.io/releases/v${nextflow_version}/nextflow-${nextflow_version}-one.jar" \
    --output "$temporary/nextflow-${nextflow_version}-one.jar"
curl "${curl_options[@]}" \
    "https://github.com/askimed/nf-test/releases/download/v${nf_test_version}/nf-test-${nf_test_version}.tar.gz" \
    --output "$temporary/nf-test.tar.gz"

printf '%s  %s\n' "$nextflow_sha256" "$temporary/nextflow" \
    "$nextflow_jar_sha256" "$temporary/nextflow-${nextflow_version}-one.jar" \
    "$nf_test_sha256" "$temporary/nf-test.tar.gz" | sha256sum --check --strict

install -m 0755 "$temporary/nextflow" "$destination/nextflow"
install -m 0644 "$temporary/nextflow-${nextflow_version}-one.jar" \
    "$nextflow_home/framework/$nextflow_version/nextflow-${nextflow_version}-one.jar"
tar -xzf "$temporary/nf-test.tar.gz" -C "$temporary" nf-test nf-test.jar
install -m 0755 "$temporary/nf-test" "$destination/nf-test"
install -m 0644 "$temporary/nf-test.jar" "$destination/nf-test.jar"
