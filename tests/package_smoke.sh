#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
smoke_root=$(mktemp -d)
trap 'rm -rf "$smoke_root"' EXIT

cd "$repo_root"
cargo package --locked --allow-dirty --target-dir "$smoke_root/target"
package_dir=$(find "$smoke_root/target/package" -mindepth 1 -maxdepth 1 -type d -name 'hap-rs-*' -print -quit)
test -n "$package_dir"
cargo test --locked --all-targets --manifest-path "$package_dir/Cargo.toml" --target-dir "$smoke_root/package-target"
cargo install --locked --path "$package_dir" --root "$smoke_root/install" --bin hap
test "$("$smoke_root/install/bin/hap" --version)" = "hap 0.1.0"
