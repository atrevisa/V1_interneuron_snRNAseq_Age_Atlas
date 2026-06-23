# V1_interneuron_snRNAseq
# V1 Interneuron snRNA-seq Age Atlas Analysis

## Overview

This repository contains the code used to generate the single-nucleus RNA sequencing (snRNA-seq) atlas of spinal V1 interneurons across postnatal development described in:

> Trevisan AJ, et al.
> *Transcriptomic analysis of spinal V1 interneurons informs their multifunctional role in motor output.*
> Nature Communications (2026)

The analysis workflow includes quality control, data integration, clustering, cell-type annotation, differential gene expression analysis, and generation of figures used in the manuscript.

---

## Associated Publication

**Citation**

Trevisan AJ, et al.

*Transcriptomic analysis of spinal V1 interneurons informs their multifunctional role in motor output.*

Nature Communications (2026)

**DOI:** To be added upon publication.

---

## Data Availability

Raw sequencing data, count matrices, and associated metadata are available through the Gene Expression Omnibus (GEO).

**GEO Accession:** GSE275595

Input data should be downloaded and placed in the appropriate data directories prior to executing the analysis scripts.

---

## Analysis Summary

The workflow performs:

1. Import and preprocessing of snRNA-seq count matrices
2. Quality-control filtering of nuclei
3. Removal of ambient RNA contamination
4. Doublet identification and exclusion
5. Normalization and scaling
6. Batch correction and dataset integration
7. Dimensionality reduction and clustering
8. Cell-type annotation and marker identification
9. Differential expression analyses
10. Generation of manuscript figures and supplementary results

---

## Repository-to-Manuscript Mapping

The table below links key analyses and manuscript figures to the scripts used to generate them.

| Manuscript Item       | Description                                                                | Script(s)                               |
| --------------------- | -------------------------------------------------------------------------- | --------------------------------------- |
| Figure 1              | Data processing, quality control, and overview of the V1 interneuron atlas | `01_preprocessing.R`                    |
| Figure 2              | Dimensionality reduction, clustering, and cell-type identification         | `02_clustering.R`                       |
| Figure 3              | Marker gene identification and subtype characterization                    | `03_marker_analysis.R`                  |
| Figure 4              | Developmental comparisons and differential expression analyses             | `04_differential_expression.R`          |
| Figure 5              | Functional characterization and downstream analyses                        | `05_downstream_analysis.R`              |
| Supplementary Figures | Additional analyses and validation experiments                             | Corresponding scripts in `scripts/`     |
| Supplementary Tables  | Marker genes and differential expression results                           | Export scripts within analysis workflow |

> **Note:** Update script names to match the final repository structure.

---

## Repository Structure

```text
.
├── data/                # Input data files (not included)
├── scripts/             # Analysis scripts
├── figures/             # Figure generation scripts and outputs
├── results/             # Processed objects and results
├── SessionInfo/         # Software versions and package information
└── README.md
```

Directory names may vary depending on the final repository version.

---

## Software Requirements

Analyses were performed using:

* R (version specified in SessionInfo)
* Seurat
* Scanpy

Additional package versions and dependency information are provided in the accompanying SessionInfo files.

---

## Running the Analysis

Scripts are intended to be run sequentially according to their filenames.

A typical workflow consists of:

1. Data preprocessing and quality control
2. Integration and clustering
3. Cell-type annotation
4. Differential expression analysis
5. Figure generation

Example:

```bash
Rscript scripts/01_preprocessing.R
Rscript scripts/02_clustering.R
Rscript scripts/03_marker_analysis.R
Rscript scripts/04_differential_expression.R
```

Refer to individual scripts for required inputs and outputs.

---

## Expected Outputs

The workflow generates:

* Processed Seurat objects
* Cluster annotations
* Marker gene tables
* Differential expression results
* Publication-quality figures
* Supplementary analysis outputs

---

## Reproducibility

All analyses were conducted using the software versions documented in the provided SessionInfo files.

Minor differences in numerical values or visualization appearance may occur when using different package versions or computing environments.

---

## Contact

Questions regarding the analysis workflow may be directed to:

**Alexandra J. Trevisan**

or

**Jay B. Bikoff**

Department of Developmental Neurobiology
St. Jude Children's Research Hospital
Memphis, TN, USA

---

## License

This repository is provided for academic and research use.

If you use this code or derived results, please cite the associated publication.
