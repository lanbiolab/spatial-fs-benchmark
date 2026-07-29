from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import torch

from spatial_fs_benchmark.integration.graphst import _build_sparse_slice_graph


def test_sparse_slice_graph_has_no_cross_slice_edges() -> None:
    coordinates = np.asarray(
        [[0, 0], [1, 0], [2, 0], [0, 0], [1, 0], [2, 0]], dtype=np.float32
    )
    slice_ids = np.asarray(["a", "a", "a", "b", "b", "b"])

    neighbors, adjacency = _build_sparse_slice_graph(coordinates, slice_ids, n_neighbors=1)

    assert neighbors.shape == (6, 6)
    assert adjacency.shape == (6, 6)
    assert neighbors.nnz == 6
    assert adjacency[0:3, 3:6].nnz == 0
    assert adjacency[3:6, 0:3].nnz == 0
    assert (adjacency != adjacency.T).nnz == 0


def test_sparse_slice_graph_handles_single_spot_slice() -> None:
    coordinates = np.asarray([[0, 0], [0, 0], [1, 0]], dtype=np.float32)
    slice_ids = np.asarray(["single", "pair", "pair"])

    neighbors, adjacency = _build_sparse_slice_graph(coordinates, slice_ids, n_neighbors=3)

    assert neighbors[0].nnz == 0
    assert adjacency[0].nnz == 0
    assert neighbors[1:3, 1:3].nnz == 2


def test_sparse_graphst_readout_matches_official_dense_operation() -> None:
    module_root = (
        Path(__file__).resolve().parents[1]
        / "external"
        / "iSTBench"
        / "Benchmark"
        / "RunModel"
        / "GraphST"
    )
    sys.path.insert(0, str(module_root))
    from GraphST.model import AvgReadout

    embedding = torch.tensor([[1.0, 2.0], [3.0, 1.0], [2.0, 4.0]])
    dense_mask = torch.tensor([[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]])
    sparse_mask = dense_mask.to_sparse().coalesce()

    readout = AvgReadout()
    dense_output = readout(embedding, dense_mask)
    sparse_output = readout(embedding, sparse_mask)

    torch.testing.assert_close(sparse_output, dense_output)
