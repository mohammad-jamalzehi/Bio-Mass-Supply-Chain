# Biomass Supply Chain Planning

**Authors:** Mohamad Jamalzehi and Kosar Ghazali.

The editable manuscript is `manuscript.tex` and its BibTeX bibliography is
`references.bib`. The compiled publication copy is
`../web/public/research/biomass-supply-chain.pdf`. The historical PDF and GAMS
source remain at the repository root.

Run `./paper/build.ps1` in PowerShell from the root, with Tectonic available.
Pass `-Tectonic <path>` if it is installed outside PATH. Tectonic downloads its
TeX bundle on the first build. The script selects the TeX Live 2022 bundle
from Tectonic's direct host. Pass `-Bundle <path-or-url>` to use a local
or alternative bundle. Tectonic runs BibTeX and the reference passes.
The script compiles into ignored `paper/build/` and copies the PDF, editable
LaTeX source and bibliography into the Vite public directory.

## Citation audit

The archived six-page PDF is dated April 21, 2023 and has no printed
bibliography or byline. It references an unnamed article. Malladi and
Sowlati's same-title 2020 paper is a verified related publication about the
UBC case. Its exact relationship to the archive has not been independently
established.
The other scholarly entries support software, solver lineage, or proposed
future methods and were checked against publisher records. The archived
report and model are listed as project artifacts. No institution or
publication venue is assigned to this manuscript.
