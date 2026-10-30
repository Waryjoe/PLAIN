# Before Publishing

- [ ] Confirm institutional and coauthor permission to share code.
- [ ] Choose the repository owner, public/private visibility and software license.
- [ ] Add verified manuscript citation and author information; no invented DOI.
- [ ] Resolve the reproducibility/provenance items in REPRODUCIBILITY.md.
- [ ] Verify final cohort selection and report exploratory subsets honestly.
- [ ] Reconcile current figure AUC points and CIs with saved computational output.
- [ ] Run full analyses in a clean environment with authorized data, record versions.
- [ ] Inspect all staged files, not only ignored/untracked files.
- [ ] Run `python3 tools/check_release.py`; check any manual additions separately.
- [ ] Check that no original data, prediction tables, model objects, case lists,
  image labels, notebook outputs, secrets, or private filesystem paths are staged.

Only initialize/version this **PLAIN-code** directory, not the enclosing study
workspace. A ZIP created before running private analyses is the clean source
snapshot; do not compress the whole directory again after adding real data.

No remote repository, commit, push or GitHub upload is performed by this package.
