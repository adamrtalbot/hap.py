#!/usr/bin/env python3
"""Validate third-party provenance, notices, and immutable fixture bytes."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = (
    Path(sys.argv[1]).resolve()
    if len(sys.argv) > 1
    else Path(__file__).resolve().parent.parent
)
MANIFEST = ROOT / "THIRD_PARTY_LICENSES" / "manifest.json"


def fail(message: str) -> None:
    raise SystemExit(f"notice policy: {message}")


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("schema") != 1:
        fail("unsupported manifest schema")

    components = manifest.get("components")
    if not isinstance(components, list) or not components:
        fail("manifest must contain at least one component")

    notices: set[str] = set()
    covered_paths: set[str] = set()
    for component in components:
        name = component.get("name", "<unnamed>")
        kind = component.get("kind")
        if kind not in {"source-port", "fixture"}:
            fail(f"{name} has unsupported component kind: {kind!r}")
        revision = component.get("revision", "")
        if len(revision) != 40 or any(
            character not in "0123456789abcdef" for character in revision
        ):
            fail(f"{name} does not use a full immutable source revision")
        source = component.get("source", "")
        if not source.startswith("https://") or revision not in source:
            fail(f"{name} source URL does not contain its immutable revision")
        if not component.get("license"):
            fail(f"{name} does not declare an SPDX license")

        notice = component.get("notice", "")
        notice_path = ROOT / notice
        if not notice_path.is_file() or not notice_path.read_text(encoding="utf-8").strip():
            fail(f"{name} notice is missing or empty: {notice}")
        notices.add(notice)

        paths = component.get("paths")
        if not isinstance(paths, list) or not paths:
            fail(f"{name} does not map any source or fixture paths")
        for relative_path in paths:
            if (
                relative_path.startswith("verification/assets/fixtures/")
                and kind != "fixture"
            ):
                fail(f"{name} fixture path must use the fixture component kind")
            if relative_path in covered_paths:
                fail(f"path is assigned more than once: {relative_path}")
            if not (ROOT / relative_path).is_file():
                fail(f"{name} source or fixture is missing: {relative_path}")
            covered_paths.add(relative_path)

        checksums = component.get("sha256", {})
        if not isinstance(checksums, dict):
            fail(f"{name} sha256 field must be an object")
        if kind == "fixture" and set(checksums) != set(paths):
            missing = sorted(set(paths) - set(checksums))
            extra = sorted(set(checksums) - set(paths))
            fail(
                f"{name} fixture checksums must cover every path "
                f"(missing={missing}, extra={extra})"
            )

        for relative_path, expected in checksums.items():
            path = ROOT / relative_path
            if relative_path not in covered_paths:
                fail(f"{name} checksum path is not listed in paths: {relative_path}")
            if not isinstance(expected, str) or len(expected) != 64 or any(
                character not in "0123456789abcdef" for character in expected
            ):
                fail(f"{name} does not use a valid SHA-256 for {relative_path}")
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            if actual != expected:
                fail(f"{name} checksum mismatch for {relative_path}")

    packaged_notices = {
        str(path.relative_to(ROOT))
        for path in (ROOT / "THIRD_PARTY_LICENSES").glob("*.txt")
    }
    if notices != packaged_notices:
        missing = sorted(packaged_notices - notices)
        stale = sorted(notices - packaged_notices)
        fail(f"manifest/notice mismatch (unmapped={missing}, missing={stale})")

    print(f"validated {len(components)} components and {len(covered_paths)} paths")


if __name__ == "__main__":
    main()
