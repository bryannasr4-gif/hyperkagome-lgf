"""
make_arxiv.py  --  build (or check) the arXiv preprint source from the journal source.

The preprint and the journal submission are the same paper.  Only the document class and the
title block differ: the journal source uses the publisher class, which prints a journal name on
page 1, and a preprint must not carry the branding of a journal that has not accepted it.

Everything from the first sectioning command to the end of the document is copied byte for byte,
and the abstract is lifted out of the journal source rather than retyped, so neither can drift.

    python make_arxiv.py            regenerate paper_arxiv/main.tex from paper/main.tex
    python make_arxiv.py --check    assert the two bodies are byte-identical; exit 1 if not

The --check mode is the drift guard.  It is falsifiable: change one character of either body,
or point it at a source built on the publisher class, and it fails.

Note on the branding test: it asks what document class the preprint declares, and whether it
references the publisher class files.  It deliberately does not search for the bare string
"ipart", which occurs inside the ordinary word "bipartite" in the body of this paper.
"""
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SRC = os.path.join(REPO, "paper", "main.tex")
DST = os.path.join(HERE, "main.tex")
SRC_FIG = os.path.join(REPO, "paper", "figs", "dos.pdf")
DST_FIG = os.path.join(HERE, "figs", "dos.pdf")

BODY_MARK = "\\section{Introduction}"
CLASS_RE = re.compile(r"\\documentclass(?:\[[^\]]*\])?\{([^}]+)\}")
PUBLISHER_MARKERS = ("ipart.cls", "ipart-layout", "{ipart}", "[cntp]")

PREAMBLE = r"""%%  arXiv preprint source.  Generated from paper/main.tex by paper_arxiv/make_arxiv.py:
%%  the preamble and title block below are the only difference; from the first sectioning
%%  command onward the source is byte-identical to the journal source.  Regenerate with
%%  "python make_arxiv.py"; verify with "python make_arxiv.py --check".
\documentclass[11pt]{article}

\usepackage{amsmath}
\usepackage{amsthm}
\usepackage{amssymb}
\usepackage{graphicx}
\usepackage{bm}
\usepackage[margin=1.1in]{geometry}
\usepackage{hyperref}

\newcommand{\rmd}{\mathrm{d}}
\newcommand{\Tr}{\mathrm{Tr}}

\theoremstyle{plain}
\newtheorem{thm}{Theorem}[section]
\newtheorem{prop}[thm]{Proposition}
\newtheorem{lem}[thm]{Lemma}
\theoremstyle{remark}
\newtheorem{rem}[thm]{Remark}

% Allow TeX a last-resort stretch rather than pushing material into the margin.
\emergencystretch=3em

\title{Lattice Green's function of the hyperkagome lattice: modular uniformization at level~30
from an orthogonal differential Galois group}

\author{%
Bryan Nasr\thanks{Independent researcher.\ \texttt{bryannasr4@gmail.com}}
\and
Jean-Marie Maillard\thanks{LPTMC, UMR 7600 CNRS, Sorbonne Universit\'e, Tour 23, 5\`eme
\'etage, case 121, 4 Place Jussieu, 75252 Paris Cedex 05, France.\
\texttt{jean-marie.maillard@sorbonne-universite.fr}}%
}

\date{}

\begin{document}

\maketitle

\begin{abstract}
%(ABSTRACT)s
\end{abstract}

\noindent\textbf{2020 Mathematics Subject Classification.}
34M15, 11F03 (Primary); 82B20, 12H05, 11F11, 11F20 (Secondary).

\medskip
\noindent\textbf{Keywords.}
Lattice Green's function; hyperkagome lattice; Picard--Fuchs operator; differential Galois
theory; symmetric square; modular curve; Fricke group; quasimodular form.

"""


def read(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return fh.read()


def split_source(text):
    """Return (abstract, body) of the journal source."""
    a = text.index("\\begin{abstract}") + len("\\begin{abstract}")
    b = text.index("\\end{abstract}")
    return text[a:b].strip("\n"), text[text.index(BODY_MARK):]


def build(text):
    abstract, body = split_source(text)
    return PREAMBLE.replace("%(ABSTRACT)s", abstract) + body


def check(src, dst):
    """Return the list of drift problems; empty means the preprint is in sync."""
    problems = []
    if BODY_MARK not in dst:
        problems.append("the preprint source has no %s" % BODY_MARK)
    elif src[src.index(BODY_MARK):] != dst[dst.index(BODY_MARK):]:
        problems.append("the two bodies differ from %s onward" % BODY_MARK)
    abstract, _ = split_source(src)
    if abstract not in dst:
        problems.append("the preprint abstract is not the journal abstract")
    cls = CLASS_RE.search(dst)
    if cls is None:
        problems.append("the preprint declares no document class")
    elif cls.group(1).strip() != "article":
        problems.append("the preprint uses the %r class, not article" % cls.group(1).strip())
    for marker in PUBLISHER_MARKERS:
        if marker in dst:
            problems.append("the preprint still references the publisher class (%r)" % marker)
    if not os.path.exists(DST_FIG):
        problems.append("figs/dos.pdf is missing, so the preprint source is not self-contained")
    elif open(DST_FIG, "rb").read() != open(SRC_FIG, "rb").read():
        problems.append("figs/dos.pdf differs from the journal figure")
    return problems


def main():
    src = read(SRC)
    if "--check" in sys.argv:
        if not os.path.exists(DST):
            print("FAIL  paper_arxiv/main.tex does not exist; run make_arxiv.py first")
            return 1
        problems = check(src, read(DST))
        if problems:
            print("FAIL  the arXiv source has drifted from the journal source:")
            for p in problems:
                print("      - %s" % p)
            return 1
        print("PASS  arXiv body byte-identical to the journal body (%d characters from %s onward)"
              % (len(src) - src.index(BODY_MARK), BODY_MARK))
        print("PASS  the abstract is the journal abstract, verbatim")
        print("PASS  document class is article; no publisher class or journal branding")
        print("PASS  figs/dos.pdf present and identical to the journal figure")
        return 0

    out = build(src)
    with open(DST, "w", encoding="utf-8", newline="") as fh:
        fh.write(out)
    os.makedirs(os.path.dirname(DST_FIG), exist_ok=True)
    shutil.copyfile(SRC_FIG, DST_FIG)
    print("wrote %s (%d characters)" % (os.path.relpath(DST, REPO), len(out)))
    print("copied %s" % os.path.relpath(DST_FIG, REPO))
    print("the arXiv package is exactly these two files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
