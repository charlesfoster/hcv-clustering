# Versioning and releases

hcv-cluster uses [semantic versioning](https://semver.org/) (`MAJOR.MINOR.PATCH`), interpreted for a
surveillance tool: the question is whether **results** can change.

| Bump | When | Examples |
|---|---|---|
| **MAJOR** | The same input and parameters can give different clusters, or an output/registry format breaks | New default threshold or region; different distance calculation; renamed or removed output columns; incompatible label-registry format |
| **MINOR** | New capability that leaves existing results unchanged | New option, output file or column, new clustering method, new subtype reference |
| **PATCH** | Bug fix or documentation that does not change results | Crash fix, clearer message, test change |

While the version is `0.x`, a MINOR bump may also contain results-affecting changes, but they are always
marked **results-affecting** in `CHANGELOG.md`. From `1.0.0` they require a MAJOR bump.

Two things need extra care because other systems depend on them:

- **Label registries are tied to a clustering definition** (metric, region, threshold, TN93 settings).
  Changing a default threshold or region therefore makes existing registries refuse to update. Treat that
  as MAJOR, say so in the changelog, and tell registry owners to start a new registry (labels from
  different definitions are not comparable).
- **Tool versions** (TN93, minimap2, MAFFT, Python) are pinned in `pixi.toml` (`core` feature). A bump
  must be checked first: run the same inputs through the old and new versions and confirm the genotypes,
  distances and clusters are identical (or record the differences as results-affecting).

## Where the version lives

`hcv_cluster/__init__.py` (`__version__`) is the source. These must match it, and
`tests/test_release.py` fails if they do not:

- `nextflow.config`: `manifest.version` and `params.pipeline_version` (the latter selects the container tag)
- `pixi.toml`: workspace `version`
- `CHANGELOG.md`: the newest `## [X.Y.Z]` heading

`pyproject.toml` reads the version from `__init__.py`.

## Making a release

1. Make sure `main` is green: `pixi run test` and `pixi run -e nf test` (the second compares the CLI with
   the Nextflow workflow).
2. Choose the version number from the table above.
3. Change the version in the four places listed, and turn the `[Unreleased]`/`unreleased` changelog entry
   into `## [X.Y.Z] - YYYY-MM-DD`, marking anything results-affecting.
4. Run the tests again, commit (`Release X.Y.Z`), and tag: `git tag -a vX.Y.Z -m "hcv-cluster X.Y.Z"`.
5. Build and publish the container with the same version as its tag:
   `docker build --build-arg GIT_COMMIT=$(git rev-parse HEAD) -t <registry>/hcv-cluster:X.Y.Z .`
   If you publish to a registry, set `process.container` (or `--pipeline_version`/a profile) accordingly.
6. Push the commit and the tag. Consumers pin the release: `nextflow run <repo> -r vX.Y.Z`.

## Reproducing a result

Every run writes `run_info.json` (pipeline version and git commit, tool versions, a hash of the genotype
panel, all parameters, and for Nextflow the engine version, profile and command line). The label registry
records the pipeline version of every run in `runs.csv`. To reproduce a result, check out the recorded tag
or commit, use the recorded container or `pixi.lock`, and rerun with the recorded parameters.
