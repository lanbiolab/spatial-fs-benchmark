# GraphST source and benchmark modifications

Original authors: Yahui Long and colleagues. Original repository and commit:
https://github.com/JinmiaoChenLab/GraphST/tree/d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74
The benchmark obtained these files through iSTBench commit `e51f60e6127dbbdd953a6fb43351b8864d9e646a`.
The pre-adaptation copies of all five Python files match the official originals.

Local sparse adaptations were recorded in the benchmark provenance on 2026-07-29:
- GraphST.py retains sparse adjacency and neighborhood matrices rather than
  converting them to dense tensors; its training input is handled accordingly.
- model.py uses sparse matrix multiplication for the neighborhood readout.
- The benchmark wrapper builds block-diagonal within-slice graphs and reduces
  reconstructed expression to 20 PCs for its common downstream evaluation.
- The original loss, optimizer and 600-epoch setting are retained.

These changes do not constitute a new GraphST algorithm. Original upstream
hashes and the local adaptation are recorded in the benchmark's
external/integration-method-references/PROVENANCE.tsv.

Revision notices were added on 2026-11-06. The explicit upstream AGPL-3.0 text is
retained as LICENSE.md; setup.py's inconsistent MIT metadata is recorded in the
audit rather than silently substituted for that text. This revision changes
source comments and licensing documentation, not the computational statements.
