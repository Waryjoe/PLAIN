"""Check public source candidates, not a substitute for institutional review."""

import argparse
import ast
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IGNORED_DIRS = {".git", "__pycache__", ".venv", "venv", "outputs", "private_workspace", "data", "models", "weights"}
TEXT_SUFFIXES = {".py", ".R", ".r", ".md", ".txt", ".yaml", ".yml"}
SPECIAL = {".gitignore", "SOURCE_MANIFEST.csv", "RELEASE_CHECKS.json", "FILE_HASHES.json"}
PATTERNS = {
    "workstation_home": re.compile(r"/(?:Users|home)/[A-Za-z0-9_.-]+/"),
    "direct_case_identifier": re.compile(r"(?<![A-Za-z0-9])(?:P|ZX)\d{6,}(?!\d)"),
    "quoted_record_number": re.compile(r"[\"']0\d{7}[\"']"),
    "labelled_record_number": re.compile(r"(?:住院号|病案号|身份证号)\s*[:：]?\s*\d{6,}"),
    "private_key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "token": re.compile(r"(?:gh[pousr]_|sk-)[A-Za-z0-9_-]{24,}"),
    "cloud_key": re.compile(r"AKIA[A-Z0-9]{16}"),
}


def public_files(root=ROOT):
    return sorted(p for p in root.rglob("*") if p.is_file()
                  and not any(x in IGNORED_DIRS for x in p.relative_to(root).parts))


def check(root=ROOT):
    errors = []
    files = public_files(root)
    for p in files:
        rel = p.relative_to(root).as_posix()
        if p.is_symlink(): errors.append(rel + ": symlink is not permitted"); continue
        if p.suffix not in TEXT_SUFFIXES and p.name not in SPECIAL:
            errors.append(rel + ": unexpected file type"); continue
        if p.stat().st_size > 2_000_000: errors.append(rel + ": unexpectedly large source file")
        try:
            text = p.read_text(encoding="utf-8-sig")
        except UnicodeError:
            errors.append(rel + ": binary or non-UTF8 file"); continue
        for name, pattern in PATTERNS.items():
            if pattern.search(text): errors.append(rel + ": " + name)
        if p.suffix == ".py":
            try: ast.parse(text)
            except SyntaxError as exc: errors.append(rel + ": " + str(exc))
    manifest = root / "SOURCE_MANIFEST.csv"
    if manifest.exists():
        with manifest.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            p = root / row["path"]
            if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest() != row["packaged_sha256"]:
                errors.append(row["path"] + ": source-manifest checksum mismatch")
    return {"public_source_files": len(files), "errors": errors,
            "privacy_scope": "Heuristic source-only scan; not clinical anonymization certification",
            "excluded_directories": sorted(IGNORED_DIRS)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result = check()
    print(json.dumps(result, indent=2) if args.json else
          f"Checked {result['public_source_files']} public source files; {len(result['errors'])} findings.\n" + "\n".join(result["errors"]))
    raise SystemExit(bool(result["errors"]))
