# E3 pilot downloads

The completed E3 implementation and exploratory pilot are published on
`codex/e3-pilot-2026-11-13`. The `main` branch is unchanged.

- [Download the complete pilot bundle (6.4 MB)](https://github.com/tarashalay88-alt/Spectral-Stability/raw/refs/heads/codex/e3-pilot-2026-11-13/deliverables/e3-pilot/E3-pilot-bundle.zip)
- [Completion and postflight report](E3-completion-report.md)
- [Scientific pilot report](E3-pilot-report.md)
- [Summary CSV](summary.csv)
- [Reproducibility manifest](run_manifest.json)
- [Implementation-only patch](e3-implementation.patch)
- [Checksums](SHA256SUMS.json)
- [Figures](figures/)

The bundle contains 100 files, including all source, 84,000 matched pilot
replication records, smoke results, tests, validation evidence, and PNG/SVG
figures. These files are also available in their experiment/source/results
locations on this branch. See `experiments/e3/README.md` for reproduction.
The frozen `experiments/e3/E3.md` is unchanged. Confirmation was not run.

The archive and historical reports were preserved byte-for-byte from the
validated local run. Statements there that no commit or remote push had been
made describe the pre-publication snapshot. Publication was subsequently
approved by the user, restricted to this new branch. The base commit at pilot
execution was `afa3dcbdb54909d82193c4a9d6f78d785d5b62db`; the implementation
source fingerprint is
`e5b01b6e772e363418487848636f9861b78bd25724304b911e7164c67c26d1c9`.
