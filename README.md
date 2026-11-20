# Cereal amylases and starch hydrolysis

Computational companion to the monograph by **Stefan Traublinger**.

[Read the German and English monograph on Zenodo](https://zenodo.org/records/24378474) · [Version DOI](https://doi.org/10.5281/zenodo.24378474)

This repository makes the selected calculations, their input files and their checks available for inspection and reuse under the component-specific conditions in [LICENSE.md](LICENSE.md). The computational baseline is research edition **1.0**, based on manuscript version **36**. Repository packaging version: **1.0.0**.

## Research question and scope

The work examines which quantities can be reconstructed from maltose-formation curves under explicitly stated kinetic and observation models. A good fit to a product curve does not, by itself, separate effective enzyme performance from operational substrate accessibility. The analytical work examines this ambiguity, the additional calibration relationships that resolve it in the specified models, and conditional parameter ranges.

The numerical examples reanalyse four digitised wheat-flour curves from two publication families. A separate INFOGEST comparison uses 20 pairs of published assay means from five laboratories. The calculations compare exponential, Weibull and reference continuations; examine sensitivity to technical assumptions; check four rational global error certificates; and test synthetic calibration and accessibility models.

The triticale section is a proposed measurement and pilot plan. These materials do not report original laboratory or baking experiments, validated transfer to triticale flour, or demonstrated improvements in bread quality or food security. The checks are internal computational checks from the same AI-assisted project, not external peer review or independent experimental validation. A historical priority claim is not established.

## Quick start

Use **Python 3.12**. From this repository directory:

```sh
python -m venv .venv
```

Activate the environment on Windows:

```bat
.venv\Scripts\activate.bat
```

Or on Linux/macOS:

```sh
source .venv/bin/activate
```

Install the recorded dependencies and run:

```sh
python -m pip install -r requirements.txt
python run_reproduction.py
```

The full run verifies the frozen research files first, then runs the calculations in a temporary copy. Results are written to `outputs/latest/`. It does not download source articles, alter the distributed scientific inputs or publish anything.

For a shorter run that rechecks stored nominal parameters:

```sh
python run_reproduction.py --skip-optimization
```

To check only the frozen baseline files and their content assignments:

```sh
python run_reproduction.py --verify-only
```

The root wrapper accepts `--output PATH` for another results directory. Do not run Python with `-O` or `-OO`, because the scientific checks use assertions.

## What is included

| Path | Purpose |
| --- | --- |
| `run_reproduction.py` | Portable entry point for the original reproduction driver. |
| `requirements.txt` | Direct dependencies recorded by the research edition. |
| `research/` | Complete frozen research edition 1.0, with its original manifests, code, inputs, tables, source register and historical document files. |
| `research/reproduce.py` | Original reproduction driver. |
| `research/verify_release.py` | Original integrity and content-assignment checks. |
| `research/zertifizierte_fehlergrenzen/` | Rational error certificates, their mathematical justification and verifier. |
| `research/DATEIRECHTE.csv` | File-specific provenance, checksums and component licensing. |
| `research/INHALTSRECHTE.csv` | Rights assignments for tables and figures in the research edition. |
| `CITATION.cff` | Citation information with the actual Zenodo monograph DOI as the preferred citation. |
| `metadata/zenodo_record.json` | Selected public record metadata verified on 19 November 2026. |
| `docs/VERIFICATION.json` | Results and limits of the new internal reproduction run. |
| `docs/PROVENANCE.md` | Relationship between the computational baseline, the book and this repository. |
| `docs/GITHUB_SCHNELLSTART_DE.md` | German instructions for uploading this prepared repository. |

The baseline retains its original directory names so that existing calculation paths and integrity checks continue to work. Files called `data/original` contain original **project input CSV files**, not original laboratory raw data or the original publications. Duplicate input copies are preserved for existing verification paths and are not additional datasets.

## Verification and limits

On 19 November 2026 the frozen baseline was rerun with nominal E/W/P2 model optimisation enabled. The delivered verification report records the executed steps, the runtime and their exit codes. The 1,572 historical technical scenarios are re-evaluated from stored parameters; they are not all independently re-optimised. Partly overlapping assertions are not independent experiments.

The run uses the available runtime with the recorded dependencies. A fresh dependency installation, GitHub Actions execution, new Word/PDF typesetting and independent human review are outside this check. See [the verification report](docs/VERIFICATION.json) and the original scientific limits in `research/README.md`.

## Publication and citation

The German and English PDFs are available at **https://zenodo.org/records/24378474**. The version DOI is **10.5281/zenodo.24378474**. The all-versions DOI is **10.5281/zenodo.24378096**.

The Zenodo record is currently titled:

> Significance and Limitations of Kinetic Models of Enzymatic Starch Hydrolysis in Cereal Systems

The preserved research edition uses the earlier title *Aussagekraft und Grenzen kinetischer Modelle der enzymatischen Stärkehydrolyse in Getreidesystemen*. `CITATION.cff` follows the currently published Zenodo record; its title was aligned on 20 November 2026. These title variants are documented explicitly rather than silently treated as byte-identical editions.

The DOI identifies the monograph record containing the two PDFs. It is **not** a DOI for this software snapshot. No software DOI, GitHub account, repository address, DNB number, ORCID or institutional affiliation has been invented. A DNB catalogue link can be added when it is available.

The original files under `research/` remain unchanged. Statements there that no DOI had yet been assigned describe the research snapshot prepared on 18 November 2026. The current publication link is documented here and in the root citation file.

## Licensing and provenance

Original project contributions are offered under **CC BY 4.0 where the relevant rights exist**. Attribution to Stefan Traublinger is required under that licence. Explicitly marked Shintani-derived components retain **CC BY-NC-SA 4.0** where protected source material is concerned. The Hackenberg entries document limited technical numerical readings and the project's own processing; they do not grant rights to the source article, figures or dissertation. INFOGEST components retain their stated CC BY 4.0 attribution requirements.

This repository is a collection with component-specific conditions, not a blanket CC BY or MIT licensing of all contents. Consult [LICENSE.md](LICENSE.md), the unchanged [research licence statement](research/LICENSE.md), and the two rights tables before reusing a particular component. The author's signature is embedded only in the monograph documents and is not separately licensed for other uses. No public-domain declaration is made.

## Deutsche Kurzfassung

Dieses Repository enthält den Rechenbestand zur Monografie von Stefan Traublinger. Die deutsche und die englische Buchfassung sind über den oben genannten Zenodo-Eintrag erreichbar. Der wissenschaftliche Ausgangsbestand unter `research/` bleibt als Forschungsausgabe 1.0 unverändert, einschließlich Quellenzuordnung, Lizenzabgrenzungen und Prüfsummen. Der neue Einstieg `python run_reproduction.py` führt die vorhandenen Prüfpfade aus. Ein Triticale-Messplan und mathematische Modellnachweise ersetzen keine eigenen Back- oder Laborversuche. DNB-Katalogdaten können später ergänzt werden.
