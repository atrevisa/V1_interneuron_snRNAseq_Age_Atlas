# V1 Interneuron snRNA-seq Age Atlas Analysis

## Overview

This repository contains the code used to generate the single-nucleus RNA sequencing (snRNA-seq) atlas of spinal V1 interneurons across postnatal development described in the preprint:

Trevisan AJ, Han K, Chapman P, Kulkarni AS, Hinton JM, Ramirez C, Klein I, Gatto G, Gabitto MI, Menon V, Bikoff JB. The transcriptomic landscape of spinal V1 interneurons reveals a role for En1 in specific elements of motor output. bioRxiv [Preprint]. 2024 Oct 26:2024.09.18.613279. doi: 10.1101/2024.09.18.613279. PMID: 39345580; PMCID: PMC11429899.

The analysis workflow includes quality control, data integration, clustering, cell-type annotation, differential gene expression analysis, and generation of figures used in the manuscript.

An interactive website can be used to view the data here: <https://v1interneurons.stjude.org/shinyApp/>


---

## Associated Publication

**Preprint:** 

Trevisan AJ, Han K, Chapman P, Kulkarni AS, Hinton JM, Ramirez C, Klein I, Gatto G, Gabitto MI, Menon V, Bikoff JB. The transcriptomic landscape of spinal V1 interneurons reveals a role for En1 in specific elements of motor output. bioRxiv [Preprint]. 2024 Oct 26:2024.09.18.613279. doi: 10.1101/2024.09.18.613279. PMID: 39345580; PMCID: PMC11429899.

**Final publication:**

In progress


---

## Data Availability

Raw sequencing data, count matrices, and associated metadata are available through the Gene Expression Omnibus (GEO).

**GEO Accession:** GSE275595

**Link to GEO page:** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275595>

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

The figures below were generated using the code in this repository

**Figure 2B-2G:** Overview of the postnatal V1 atlas, clustering at multiple resolutions, identification of previously characterized non-overlapping clades (FoxP2+, Pou6f2+, Sp8+, and Calb1+)

**Figure 3A, 3H, 3I:** Identification of a novel V1 clade identified by Rnf220 expression that is mutually exclusive to FoxP2, Pou6f2, Sp8, and Calb1

**Figure 4A, 4C, 4E:** More detailed molecular profiling of V1 interneurons, specifically, the differentially expressed transcription factors and ion channels

**Figure 5A-5G:** Analysis of age-related changes in V1 gene expression

**Supplementary Figure 2D - 2N:** Detailed description of how V1 nuclei were identified including QC, positive, and negative selection

**Supplementary Figure 3:** Basic QC stats of the final data included in the figures and website

**Supplementary Figure 5A - 5D:** Detailed description of previously characterized transcription factors within the V1 interneurons, and detailed profiling of the Sp8+ subsets

**Supplementary Figure 6:** Analyses showing the relationships between the V1 subsets

**Supplementary Table 1:** QC metrics from cellranger

**Supplementary Table 2:** Marker genes used in V1 nuclei identification

**Supplementary Table 3:** The unique cell barcodes that were used in the final V1 atlas

**Supplementary Table 4:** Differentially expressed genes in V1 subsets

**Supplementary Table 6:** Differentially expressed genes across time points


---

## Software Requirements

Analyses were performed using:

* R (version specified in SessionInfo)
* Python
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

Refer to individual scripts for required inputs and outputs.

---

## Expected Outputs

The workflow generates:

* Processed Seurat objects
* Cluster annotations
* Marker gene tables
* Differential expression results


---

## Reproducibility

All analyses were conducted using the software versions documented in the provided SessionInfo files.

Minor differences in numerical values or visualization appearance may occur when using different package versions or computing environments.

---

## Contact

Questions regarding the analysis workflow may be directed to:

**Jay B. Bikoff**

<Jay.Bikoff@STJUDE.ORG>  
Department of Developmental Neurobiology  
St. Jude Children's Research Hospital  
Memphis, TN, USA

---

## License

This repository is provided for academic and research use.

If you use this code or derived results, please cite the associated publication.
