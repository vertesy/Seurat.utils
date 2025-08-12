# AGENTS

## Overview
Seurat.utils is an R package with helper functions that extend [Seurat](https://satijalab.org/seurat) for single-cell analysis. It collects routines for metadata handling, visualization, and general utilities.

## Repository layout
- `R/` – package source code.
  - `Seurat.Utils.R` – core functions.
  - `Seurat.Utils.Metadata.R` – metadata helpers.
  - `Seurat.Utils.Visualization.R` – plotting helpers.
  - `Seurat.utils.less.used.R` – rarely used functions.
  - `*.bac` files are backups; ignore them.
- `man/` – function documentation generated from roxygen comments.
- `Vignettes/` – long-form usage examples.
- `Example.usage/` – short code examples.
- `Development/` – experimental scripts.
- `Function.Dependencies.md` – overview of internal function relationships.

## Working with the code
1. Place new or modified functions in `R/` and document them with roxygen2 comments.
2. Run package checks before committing:
   ```bash
   R -q -e "devtools::document(); devtools::check()"
   ```
   This updates documentation and performs `R CMD check`.
3. Commit only `.R` and generated documentation files; do not commit `.bac` backups or temporary files.

## Coding style
- Prefer R's native pipe `|>` rather than `%>%`.
- Follow the tidyverse conventions and the styles used in other @vertesy packages.

## Dependencies
Core functionality expects several other @vertesy packages to be installed:
- [Stringendo](https://github.com/vertesy/Stringendo)
- [CodeAndRoll2](https://github.com/vertesy/CodeAndRoll2)
- [ReadWriter](https://github.com/vertesy/ReadWriter)
- [MarkdownHelpers](https://github.com/vertesy/MarkdownHelpers)
- [MarkdownReports](https://github.com/vertesy/MarkdownReports)
- [ggExpress](https://github.com/vertesy/ggExpress)
Optional but recommended:
- [DatabaseLinke.R](https://github.com/vertesy/DatabaseLinke.R)
- [Rocinante](https://github.com/vertesy/Rocinante)

Ensure these packages are installed when running or testing the code.

## Getting started
New contributors should read `README.md` for installation instructions and explore `R/Seurat.Utils.R` to see typical function style. `Function.Dependencies.md` and the vignettes provide deeper insights into the package's structure.

Happy hacking!
