#!/usr/bin/env python3
"""Require zlib-rs to be the only resolved DEFLATE implementation."""

from __future__ import annotations

import json
import subprocess


APPROVED_BACKEND = "zlib-rs"
NATIVE_TARGETS = (
    "aarch64-apple-darwin",
    "aarch64-unknown-linux-gnu",
    "x86_64-apple-darwin",
    "x86_64-pc-windows-msvc",
    "x86_64-unknown-linux-gnu",
    "x86_64-unknown-linux-musl",
)
BANNED_BACKENDS = {
    "cloudflare-zlib-sys",
    "libdeflate-sys",
    "libdeflater",
    "libz-ng-sys",
    "libz-sys",
    "miniz_oxide",
}


def main() -> None:
    resolved_names: set[str] = set()
    for target in NATIVE_TARGETS:
        result = subprocess.run(
            [
                "cargo",
                "metadata",
                "--locked",
                "--format-version",
                "1",
                "--filter-platform",
                target,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        metadata = json.loads(result.stdout)
        resolved_ids = {node["id"] for node in metadata["resolve"]["nodes"]}
        resolved_names.update(
            package["name"]
            for package in metadata["packages"]
            if package["id"] in resolved_ids
        )

    active_banned = sorted(resolved_names & BANNED_BACKENDS)
    if active_banned:
        raise SystemExit(
            "compression policy: alternate DEFLATE backend(s) resolved: "
            + ", ".join(active_banned)
        )
    if APPROVED_BACKEND not in resolved_names:
        raise SystemExit(
            f"compression policy: approved backend {APPROVED_BACKEND} is not resolved"
        )

    print(f"compression policy: {APPROVED_BACKEND} is the sole resolved backend")


if __name__ == "__main__":
    main()
