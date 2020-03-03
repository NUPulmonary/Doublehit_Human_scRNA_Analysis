# An analysis workflow for integrating human lung scRNA-seq data to investigate age-related heterogeneity
## Abstract

(TEMP) A dysfunctional response to inhaled pathogens and toxins drives a substantial portion of the susceptibility to acute and chronic lung disease in the elderly. Using genetic lineage tracing, heterochronic adoptive transfer, parabiosis, and treatment with metformin, we found the lung microenvironment drives age-related transcriptomic changes in alveolar macrophages that include reductions in cell cycle genes and increased expression of inflammatory genes.  These changes are independent of alveolar macrophage ontogeny, circulating factors or circulating monocytes.  Changes in the microenvironment, including changes in extracellular matrix composition, induce a resistance to proliferative signals from CSF2. Severe injury can induce the replacement of long-lived tissue resident alveolar macrophages with monocyte-derived alveolar macrophages, but both respond similarly to a subsequent injury.  These findings place the lung microenvironment upstream of the dysfunctional immune responses to inhaled environmental challenge in aging.  

Single cell RNA-seq (scRNA-seq) captures the transcriptomic phenotype of multiple cell populations within a tissue simultaneously. We utilized widely used R package “Seurat” and Canonical Correlation Analysis procedure to aggregate and analyze together data from 6 published dataset. Our integration included a total number of 38 samples, covering age from 17 to 88. The merged dataset provided sufficient statistical power and homogeneity to allow discovery of common aging biomarkers across distinct cell populations. We concluded that there were no heterogeneity or emerging new cell groups in avalor macrophages, which is consistent with our observation in mouse. Through pseudo bulk, we identified 673 differentially expressed genes between young and old samples and these genes were significantly overlapped with our bulk RNAseq analysis in mouse alveolar macrophages.

## Prerequisites

Install all required R packages in the R_requirement.txt files using either bioconductor or CRAN

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(R_requirement.txt)
```
or

```
install.packages(R_requirement.txt)
```

## Data availability
We included 6 public available scRNA-seq datasets from lungs of healthy controlled or donor in our analysis workflow. The source of data are list below:

* **Reyfman et al. 2019**  - [Single-Cell Transcriptomic Analysis of Human Lung Provides Insights into the Pathobiology of Pulmonary Fibrosis](https://www.ncbi.nlm.nih.gov/pubmed/30554520)
* **Madissoon et al. 2020** - [scRNA-seq assessment of the human lung, spleen, and esophagus tissue stability after cold preservation](https://www.ncbi.nlm.nih.gov/pubmed/31892341)
* **Raredon et al. 2019** - [Single-cell connectomic analysis of adult mammalian lungs](https://www.ncbi.nlm.nih.gov/pubmed/31840053)
* **Morse et al. 2019** -[Proliferating SPP1/MERTK-expressing macrophages in idiopathic pulmonary fibrosis](https://www.ncbi.nlm.nih.gov/pubmed/31221805)
* **Habermann et al. 2019** -[Single-cell RNA-sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis](https://www.biorxiv.org/content/10.1101/753806v1)
* **Valenzi et al. 2019** -[Single-cell analysis reveals fibroblast heterogeneity and myofibroblasts in systemic sclerosis-associated interstitial lung disease](https://www.ncbi.nlm.nih.gov/pubmed/31405848)


## Results

This workflow includes a R markdown file to guide the readers step by step for our analysis workflow. The following only highlights some of the key findings:

###The general workflow of our analysis is:
![] (resources/flowchart1.png)

###The total number samples is 52 and after QC control, the number of samples after filtering is 38
![] (resources/samplefiltering.png)

### We use standard Seurat [SCTtransform](https://satijalab.org/seurat/v3.1/integration.html) pipeline to perform integration. After integration, there are a total number of 

```

```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Versioning

We use Seurat V3.1.2 and R V3.5.1. on Northwestern High Performance Computing Cluster

## Authors

* **Ziyou Ren** - *Phd student in Bioinformatics* - [Northwestern University](https://amaral.northwestern.edu/people/ren/)
* **More authors**

See also the list of [contributors](https://github.com/NUPulmonary/Doublehit_Human_scRNA_Analysis/commits) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* We thank Northwestern IT and QUEST for their support
* NIH funding
* Driskill Graduate Program in Life Science