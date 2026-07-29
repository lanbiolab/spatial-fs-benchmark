from spatial_fs_benchmark.feature_selection.hybrid import HybridHVGSVGSelector


def test_balanced_union_fills_overlap_to_requested_size() -> None:
    selected = HybridHVGSVGSelector._balanced_union(
        ["a", "b", "c", "d", "e", "f"],
        ["a", "b", "x", "y", "z", "q"],
        target=6,
    )
    assert len(selected) == 6
    assert len(set(selected)) == 6
    assert {"a", "b", "c", "x"}.issubset(selected)


def test_intersection_is_exact_and_ranked_by_both_lists() -> None:
    selected = HybridHVGSVGSelector._intersection(
        ["a", "b", "c", "d"],
        ["c", "b", "x", "a"],
        target=4,
    )
    assert selected == ["b", "c", "a"]
