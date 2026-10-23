import csv
from pathlib import Path

import pytest

from hcv_cluster import labels

DEF = {"metric": "tn93", "region": "e1-e2", "threshold": 0.03}


def _partition(*groups: tuple[str, ...], genotype: str = "1a") -> dict[tuple[str, str], str]:
    """Each group of sample names is one run-local component."""
    return {
        (genotype, sample): f"{genotype}_C{index:04d}"
        for index, group in enumerate(groups, start=1)
        for sample in group
    }


def _run(registry, partition, date="2026-08-05", run_id="r", track="tn93", definition=DEF):
    return labels.assign_labels(
        registry, track, definition, partition, run_id=run_id, run_date=date
    )


def test_new_clusters_are_named_by_month_and_singletons_are_not() -> None:
    registry = labels.Registry()
    # A compact ISO date must be normalised, not sliced as if it were YYYY-MM-DD.
    result = _run(registry, _partition(("a", "b"), ("c",), ("d", "e", "f")), date="20260805")
    assert registry.labels["HCV-2608-001"].created_date == "2026-08-05"
    assert result.labels[("1a", "a")] == result.labels[("1a", "b")] == "HCV-2608-001"
    assert result.labels[("1a", "d")] == "HCV-2608-002"
    assert result.labels[("1a", "c")] == ""
    assert result.new_labels == {"HCV-2608-001", "HCV-2608-002"}


def test_counter_restarts_each_month_but_never_reuses_a_number() -> None:
    registry = labels.Registry()
    _run(registry, _partition(("a", "b")), date="2026-08-05")
    result = _run(registry, _partition(("a", "b"), ("x", "y")), date="2026-08-20", run_id="r2")
    assert result.labels[("1a", "x")] == "HCV-2608-002"
    result = _run(registry, _partition(("a", "b"), ("x", "y"), ("p", "q")), date="2026-09-01", run_id="r3")
    assert result.labels[("1a", "p")] == "HCV-2609-001"


def test_label_survives_growth_and_renumbered_run_local_ids() -> None:
    registry = labels.Registry()
    _run(registry, _partition(("a", "b")))
    # "0" sorts first, so the run-local ID of the existing cluster changes; its label must not.
    result = _run(registry, _partition(("0", "z"), ("a", "b", "c")), run_id="r2")
    assert result.labels[("1a", "a")] == result.labels[("1a", "c")] == "HCV-2608-001"
    assert result.labels[("1a", "0")] == "HCV-2608-002"


def test_merge_keeps_the_oldest_label_and_records_the_others() -> None:
    registry = labels.Registry()
    _run(registry, _partition(("a", "b"), ("c", "d")))
    result = _run(registry, _partition(("a", "b", "c", "d", "e")), run_id="r2")
    assert set(result.labels.values()) == {"HCV-2608-001"}
    assert registry.labels["HCV-2608-002"].status == "merged"
    assert registry.labels["HCV-2608-002"].merged_into == "HCV-2608-001"
    assert registry.membership[("tn93", "1a", "c")] == "HCV-2608-001"
    assert [e["event"] for e in result.events if e["event"] == "merged"] == ["merged"]
    assert registry.resolve("HCV-2608-002") == "HCV-2608-001"


def test_split_keeps_label_with_largest_piece_and_names_the_rest() -> None:
    registry = labels.Registry()
    _run(registry, _partition(("a", "b", "c", "d", "e")))
    result = _run(registry, _partition(("a", "b", "c"), ("d", "e")), run_id="r2", date="2026-09-02")
    assert result.labels[("1a", "a")] == "HCV-2608-001"
    assert result.labels[("1a", "d")] == "HCV-2609-001"
    assert registry.labels["HCV-2609-001"].split_from == "HCV-2608-001"
    assert registry.membership[("tn93", "1a", "d")] == "HCV-2609-001"


def test_absent_members_do_not_lose_the_label_and_a_lone_survivor_is_unlabelled() -> None:
    registry = labels.Registry()
    _run(registry, _partition(("a", "b", "c")))
    result = _run(registry, _partition(("a",), ("x", "y")), run_id="r2")
    assert result.labels[("1a", "a")] == ""
    assert registry.membership[("tn93", "1a", "a")] == "HCV-2608-001"
    assert registry.membership[("tn93", "1a", "b")] == "HCV-2608-001"
    # a returns to a cluster later: it brings its label back
    result = _run(registry, _partition(("a", "b")), run_id="r3")
    assert result.labels[("1a", "a")] == "HCV-2608-001"


def test_genotypes_are_independent_and_numbers_are_shared_across_tracks() -> None:
    registry = labels.Registry()
    part = {**_partition(("a", "b"), genotype="1a"), **_partition(("a", "b"), genotype="3a")}
    result = _run(registry, part)
    assert result.labels[("1a", "a")] != result.labels[("3a", "a")]
    other = _run(
        registry, _partition(("a", "b")), track="timed", definition={"metric": "timed"}, run_id="r2"
    )
    assert other.labels[("1a", "a")] == "HCV-2608-003"
    assert registry.membership[("tn93", "1a", "a")] == "HCV-2608-001"


def test_definition_change_is_refused() -> None:
    registry = labels.Registry()
    _run(registry, _partition(("a", "b")))
    with pytest.raises(labels.RegistryError, match="different registry"):
        _run(registry, _partition(("a", "b")), definition={**DEF, "threshold": 0.045})


def test_label_run_persists_only_when_committed(tmp_path: Path) -> None:
    kwargs = dict(run_id="r1", run_date="2026-08-05", outdir="out")
    registry_dir = tmp_path / "registry"
    result, table, _ = labels.label_run(
        registry_dir, "tn93", DEF, _partition(("a", "b")), commit=False, **kwargs
    )
    assert result.labels[("1a", "a")] == labels.UNASSIGNED
    assert {e["label"] for e in result.events} == {labels.UNASSIGNED}
    assert not registry_dir.exists() or not (registry_dir / labels.LABELS_FILE).exists()

    result, table, _ = labels.label_run(
        registry_dir, "tn93", DEF, _partition(("a", "b")), commit=True, **kwargs
    )
    assert result.labels[("1a", "a")] == "HCV-2608-001"
    assert table[0]["n_members_total"] == 2

    reloaded = labels.load_registry(registry_dir)
    assert reloaded.membership[("tn93", "1a", "a")] == "HCV-2608-001"
    assert reloaded.tracks["tn93"]
    assert [r["event"] for r in reloaded.events] == ["created"]

    # A read-only run now sees the saved label and mints nothing.
    result, _table, _ = labels.label_run(
        registry_dir, "tn93", DEF, _partition(("a", "b", "c")), commit=False,
        run_id="r2", run_date="2026-09-01", outdir="out",
    )
    assert result.labels[("1a", "c")] == "HCV-2608-001"
    assert len(list(csv.DictReader((registry_dir / labels.RUNS_FILE).open()))) == 1


def test_split_provenance_does_not_depend_on_sample_name_order() -> None:
    for names in (("a4", "a5", "a6", "a7"), ("r4", "r5", "r6", "r7")):
        registry = labels.Registry()
        _run(registry, _partition(("q1", "q2", "q3", "q4")))
        _run(registry, _partition(("q1", "q2", "q3", "q4"), ("p1", "p2", "p3", *names)), run_id="r2")
        assert registry.labels["HCV-2608-002"].created_seq == 2
        result = _run(
            registry,
            _partition(("p1", "p2", "p3", "q1", "q2", "q3", "q4"), names[:2], names[2:]),
            run_id="r3",
        )
        splits = [e for e in result.events if e["event"] == "split"]
        assert {e["related_label"] for e in splits} == {"HCV-2608-002"}


def test_absorbed_labels_follow_merge_chains() -> None:
    registry = labels.Registry()
    _run(registry, _partition(("z1", "z2"), ("a1", "a2"), ("b1", "b2")))
    _run(registry, _partition(("z1", "z2"), ("a1", "a2", "b1", "b2")), run_id="r2")
    result = _run(registry, _partition(("z1", "z2", "a1", "a2", "b1", "b2")), run_id="r3")
    row = labels.label_table(registry, "tn93", result.labels)[0]
    assert row["label"] == "HCV-2608-001"
    assert row["absorbed_labels"] == "HCV-2608-002;HCV-2608-003"


def test_plots_name_clusters_by_persistent_label_and_keep_run_local_ids_otherwise() -> None:
    from hcv_cluster import viz

    merged = [
        {"sample_id": "a", "cluster_id": "1a_C0001", "cluster_label": "HCV-2608-001", "cluster_size": "2"},
        {"sample_id": "n", "cluster_id": "1a_C0002", "cluster_label": labels.UNASSIGNED, "cluster_size": "2"},
        {"sample_id": "s", "cluster_id": "1a_C0003", "cluster_label": "", "cluster_size": "1"},
    ]
    shown = {r["sample_id"]: r["cluster_id"] for r in viz.apply_persistent_labels(merged)}
    assert shown == {"a": "HCV-2608-001", "n": "1a_C0002", "s": "1a_C0003"}

    # A per-genotype table has no label column and borrows labels from the merged one.
    per_genotype = [{"sample_id": "a", "cluster_id": "C0001", "cluster_size": "2"}]
    shown = viz.apply_persistent_labels(per_genotype, viz.persistent_label_map(merged))
    assert shown[0]["cluster_id"] == "HCV-2608-001"
