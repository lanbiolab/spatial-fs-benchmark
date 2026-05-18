# Metric Alignment Notes

This note documents how the current benchmark aligns with the metric system used
in `iSTBench`.

## Metrics implemented

The current result table includes the following paper-style metrics:

- integration: `dASW`, `dLISI`, `ILL`, `bASW`, `iLISI`, `GC`
- clustering: `ari`, `nmi`, `CHAOS`, `PAS`
- alignment: `Accuracy`, `Ratio`
- slice representation: `slice_repr_ARI`, `slice_repr_NMI`

## Alignment status

### Directly aligned in spirit and naming

- `CHAOS`, `PAS`
- `Accuracy`, `Ratio`
- abundance-based slice representation `ARI/NMI`

### Close approximations of the iSTBench/scIB stack

The paper computes integration metrics through `scib_metrics.Benchmarker`. The
current code reproduces the same metric family and naming, but uses local
implementations instead of the external `scib_metrics` package:

- `dASW`
- `dLISI`
- `ILL`
- `bASW`
- `iLISI`
- `GC`

This keeps the benchmark self-contained and easier to maintain, but the values
should be treated as benchmark-compatible approximations rather than byte-for-
byte reproductions of the official `scib_metrics` outputs.

## Current dataset limitation

The sample `BaristaSeq` data used for the MVP run does not contain a
`slice_class` ground-truth column. Because of that:

- `slice_repr_ARI` and `slice_repr_NMI` are present in the framework
- on `BaristaSeq`, they may be missing or uninformative unless slice-level class
  labels are added in dataset metadata

## Recommended next step for stricter paper replication

If exact paper-level parity is required, the next extension should be:

1. add an optional adapter around `scib` / `scib_metrics`
2. allow dataset-specific `slice_class_key`
3. support domain-count sweeps for slice representation, mirroring the original
   `4..10` domain-count evaluation
