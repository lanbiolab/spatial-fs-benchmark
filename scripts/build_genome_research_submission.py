#!/usr/bin/env python3

"""Build a clean Genome Research upload package from the current manuscript."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def run(command: list[str], cwd: Path | None = None) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--manuscript-dir",
        type=Path,
        default=Path("manuscript_genome_research_spatial_omics_v2"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("Genome_Research_Submission_Upload_Files"),
    )
    args = parser.parse_args()

    root = args.project_root.resolve()
    manuscript = (root / args.manuscript_dir).resolve()
    output = (root / args.output_dir).resolve()
    pandoc = root / ".conda-env" / "bin" / "pandoc"
    reference_doc = root / "manuscript_genome_research_special_issue" / "genome_research_reference.docx"

    required = [
        manuscript / "main.tex",
        manuscript / "main.pdf",
        manuscript / "supplement.tex",
        manuscript / "supplement.pdf",
        manuscript / "references.bib",
        pandoc,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs: " + ", ".join(missing))

    output.mkdir(parents=True, exist_ok=True)
    figures_out = output / "optional_separate_figures" / "main"
    figures_out.mkdir(parents=True, exist_ok=True)

    # Current PDF manuscript and standalone main figures.
    shutil.copy2(manuscript / "main.pdf", output / "01_manuscript_Genome_Research.pdf")
    for index in range(1, 7):
        source = manuscript / "figures" / f"Figure{index}.pdf"
        if not source.exists():
            raise FileNotFoundError(source)
        shutil.copy2(source, figures_out / f"Figure{index}.pdf")

    # Pandoc embeds the PDF figures, but PNG companions make Word rendering
    # consistent across desktop and browser-based Word clients.
    word_figures = manuscript / "word_figures_current"
    word_figures.mkdir(exist_ok=True)
    for index in range(1, 7):
        run(
            [
                "pdftoppm",
                "-singlefile",
                "-png",
                "-r",
                "240",
                str(manuscript / "figures" / f"Figure{index}.pdf"),
                str(word_figures / f"Figure{index}"),
            ]
        )

    tex = (manuscript / "main.tex").read_text(encoding="utf-8")
    frontmatter = r"""\begin{document}
\begin{center}
\textbf{Feature selection objectives shape spatial signal recovery and downstream analysis in spatial transcriptomics}

Wei Lan\textsuperscript{1,*}, Yang Yan\textsuperscript{1}, Xuhua Yan\textsuperscript{1}, Qingfeng Chen\textsuperscript{1}, Min Li\textsuperscript{2}, and Yi Pan\textsuperscript{3}
\end{center}

\noindent\textsuperscript{1}School of Computer, Electronics and Information, Guangxi University, Nanning 530004, Guangxi, China.

\noindent\textsuperscript{2}School of Computer Science and Engineering, Central South University, Changsha 410083, Hunan, China.

\noindent\textsuperscript{3}Shenzhen University of Advanced Technology, Shenzhen 518107, Guangdong, China.

\noindent\textsuperscript{*}Correspondence: Wei Lan; lanwei@gxu.edu.cn
"""
    marker = "\\begin{document}\n\\maketitle"
    if marker not in tex:
        raise RuntimeError("Could not locate the manuscript title block")
    tex = tex.replace(marker, frontmatter, 1).replace("\\linenumbers", "")
    # Pandoc treats the LaTeX abstract environment as document metadata and
    # otherwise moves it ahead of the explicit title/author block.
    tex = tex.replace("\\begin{abstract}", "\\section*{Abstract}", 1)
    tex = tex.replace("\\end{abstract}", "", 1)
    for index in range(1, 7):
        tex = tex.replace(
            f"{{Figure{index}.pdf}}",
            f"{{word_figures_current/Figure{index}.png}}",
        )
    docx_source = manuscript / "docx_source_current.tex"
    docx_source.write_text(tex, encoding="utf-8")

    pandoc_command = [
        str(pandoc),
        str(docx_source),
        "--from=latex",
        "--to=docx",
        "--citeproc",
        f"--bibliography={manuscript / 'references.bib'}",
        f"--resource-path={manuscript}:{word_figures}",
        f"--reference-doc={reference_doc}",
        "-o",
        str(output / "01_manuscript_Genome_Research.docx"),
    ]
    run(pandoc_command)

    # The portal accepts pasted cover-letter text, but a current editable copy
    # is retained in the package for audit and reuse.
    run(
        [
            str(pandoc),
            str(manuscript / "cover_letter.tex"),
            "--from=latex",
            "--to=docx",
            f"--reference-doc={reference_doc}",
            "-o",
            str(output / "02_cover_letter_Genome_Research.docx"),
        ]
    )

    shutil.copy2(manuscript / "supplement.pdf", output / "04_supplemental_material.pdf")

    source_data_zip = output / "05_supplemental_source_data.zip"
    source_data_zip.unlink(missing_ok=True)
    run(
        ["zip", "-q", "-r", str(source_data_zip), "source_data"],
        cwd=manuscript,
    )

    latex_zip = output / "03_Genome_Research_LaTeX_source.zip"
    latex_zip.unlink(missing_ok=True)
    zip_inputs = [
        "main.tex",
        "supplement.tex",
        "references.bib",
        "authblk.sty",
        "figures",
        "supplemental_figures",
        "supplemental_tables",
    ]
    run(["zip", "-q", "-r", str(latex_zip), *zip_inputs], cwd=manuscript)


if __name__ == "__main__":
    main()
