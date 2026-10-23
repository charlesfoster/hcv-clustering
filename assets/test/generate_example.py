#!/usr/bin/env python3
"""Regenerate assets/test/example.fasta: a small synthetic HCV 1a data set.

Everything is derived from the public subtype 1a reference record in
reference_data/genbank (no patient data), with fixed random seeds. It contains a
5-sequence cluster, a pair, singletons, a mixed infection, an exact repeat, a
sequence with too much missing data and a non-HCV sequence, so one run exercises
every exclusion stage of the pipeline.
"""
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import hcv_cluster_prep as prep  # noqa: E402

spec = prep.DEFAULT_REFERENCE_CATALOG["1a"]
reference = prep.read_fasta_records(prep.DEFAULT_REFS_DIR / f"{spec.accession}.fasta")[0].sequence
region = reference[900:3100]  # envelope genes plus flanks, ~2.2 kb

rng = random.Random(7)


def mutate(sequence: str, n: int, rng: random.Random) -> str:
    chars = list(sequence)
    for position in rng.sample(range(len(chars)), n):
        chars[position] = rng.choice([b for b in "ACGT" if b != chars[position]])
    return "".join(chars)


founder = mutate(region, 40, rng)          # unrelated to the other lineages below
lineage_b = mutate(region, 90, rng)
records: list[tuple[str, str]] = []
for i in range(5):
    records.append((f"clusterA_{i + 1}", mutate(founder, 2, rng)))
for i in range(2):
    records.append((f"pairB_{i + 1}", mutate(lineage_b, 2, rng)))
records.append(("single_1", mutate(region, 160, rng)))
records.append(("single_2", mutate(region, 200, rng)))
records.append(("mixed_1", mutate(founder, 3, rng)))                # mixed infection:
records.append(("mixed_1", mutate(lineage_b, 3, rng)))              #   two distinct sequences
records.append(("clusterA_1", records[0][1]))                       # exact repeat
low = mutate(founder, 2, rng)
records.append(("low_coverage", low[: len(low) // 2] + "N" * (len(low) - len(low) // 2)))
records.append(("not_hcv", "".join(rng.choice("ACGT") for _ in range(1500))))

out = Path(__file__).with_name("example.fasta")
out.write_text("".join(f">{h}\n{s}\n" for h, s in records))
print(f"wrote {len(records)} records to {out}")
