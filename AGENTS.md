# Guidance for Agents

Version: 2026.08.31-00:00

## I: Generic (all @vertesy repos)

### 3. General rules 

1. Write condense and very cleary understandable code and text.
2. ALWAYS understand the larger goal and the full context first: read the whole file, grep all call sites, check the roxygen block.
3. Use short inline comments to explain the code: separate lines before a larger block for the higher-level intent, and/or short trailing comments after a line for that specific step.

- Every function starts with a COMPACT input-argument assertion for key inputs. Use combined stopifnot statements.

### 2. Code Review Rules

Make every finding easy to scan and understand.

- Use simple, direct English.
- Use short sentences and bullet points.
- Avoid compiler jargon, dense technical language, and noun stacking.
- Name files, functions, variables, and arguments explicitly.
- Never use vague references such as "the new formal" or "subsequent values".

For each finding, state:

- **Problem:** What is wrong or will break?
- **Trigger:** When does it happen?
- **Fix:** What should be changed?

Keep only findings that can be explained clearly and concisely.
Do not flag formatting, line length, or missing tests.

### 3. Update the Source, Not Just the Documentation

Documentation is generated from upstream sources: `.Rd` files from roxygen annotations and `DESCRIPTION` from `Development/Dependencies.R` via `config.R`.

Package rebuilds overwrite these files, so always update the upstream source first, then regenerate the documentation.

## II: Repos of R function libraries

- New arguments go at the end, just before `...`. Never insert in the middle.
- Do not use tests.
- Never update the package version unless the user explicitly requests a version change.
- Do not raise code review findings that ask for a package version change.

## III: Seurat.utils specific

**Seurat.utils** — R helper functions extending Seurat for single-cell analysis: metadata handling, visualization, general utilities.

- `R/`: `Seurat.Utils.R` (core), `Seurat.Utils.Metadata.R` (metadata), `Seurat.Utils.Visualization.R` (plotting), `Seurat.utils.less.used.R` (rarely used). `*.bac` files are backups — ignore, never commit.
- `Function.Dependencies.md`: map of internal function relationships.
- Prefer native pipe `|>` over `%>%`.
- Depends on @vertesy `Stringendo`, `CodeAndRoll2`, `ReadWriter`, `MarkdownHelpers`, `MarkdownReports`, `ggExpress`. Optional: `DatabaseLinke.R`, `Rocinante`.
