import ast
import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from plain_paths import plain_path


class PackageTests(unittest.TestCase):
    def test_all_python_parses(self):
        for p in ROOT.rglob("*.py"):
            if "private_workspace" not in p.parts:
                ast.parse(p.read_text())

    def test_default_paths_are_local(self):
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(plain_path("PROJECT"), str(ROOT / "private_workspace/project"))
            self.assertEqual(plain_path("REPO", "scripts"), str(ROOT / "scripts"))

    def test_override_and_trailing_slash(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, {"PLAIN_PROJECT_ROOT": tmp}):
            self.assertEqual(plain_path("PROJECT", "abc/"), str(Path(tmp) / "abc") + os.sep)

    def test_unknown_root_fails(self):
        with self.assertRaises(ValueError): plain_path("UNKNOWN")

    def test_source_manifest_and_privacy(self):
        spec = importlib.util.spec_from_file_location("checker", ROOT / "tools/check_release.py")
        checker = importlib.util.module_from_spec(spec); spec.loader.exec_module(checker)
        self.assertEqual(checker.check()["errors"], [])

    def test_data_files_ignored(self):
        patterns = (ROOT / ".gitignore").read_text().splitlines()
        for pattern in ("private_workspace/", "outputs/", "*.xlsx", "*.rds", "*.nii.gz", "*.svs", "*.ipynb"):
            self.assertIn(pattern, patterns)


if __name__ == "__main__":
    unittest.main()
