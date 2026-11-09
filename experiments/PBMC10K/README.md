## PBMC10K Data

### 1. Data preparation

We used the processed, cell-type-annotated PBMC10K dataset from [TIGER (Zenodo)](https://zenodo.org/records/7425777). Download the data and rename the extracted directory to `TIGER_zenodo`.

**Input data:** `TIGER_zenodo/data/PBMC/pbmc_final.rds`

We used [CollecTRI](https://github.com/saezlab/CollecTRI) as the prior gene regulatory network (GRN) for all cell types.

Generate the necessary input files:

```bash
cp PBMC10K_v2.R TIGER_zenodo/code
Rscript TIGER_zenodo/code/PBMC10K_v2.R
```

**Conda environment:**

```text
The conda environment for running TIGER, plus additional R dependencies.
```

The script generates cell names, cell-type-specific expression matrices, regulator and gene lists, and TFAct input files. Search for `write.table` in the script for the complete list of outputs.

### 2. Run TFA inference algorithms

Execute the Snakemake workflow (ensure that the matlab modules, paths to conda environment in `rules/<algorithm>.smk`, and config.yml are correctly specified as in the simulation README.md):

```bash
snakemake -k --cores 50 --use-conda --rerun-incomplete
```

The workflow produces outputs similar to those in the simulation experiments, but no ground-truth evaluation is performed.

### 3. Compare TFA with chromVAR

After Snakemake completes, run:

```bash
Rscript get_correlation.R
```

This computes Pearson and Spearman correlations between inferred TFA and chromVAR deviation scores.

**Outputs:**

| Analysis | Output |
|---|---|
| Per-TF correlation summary across all cells (mean and SD) | `TFA_output/analysis/TFA_corr/overall/mean_(Spearman\|Pearson).pdf` |
| Correlation coefficient distributions (boxplots) | `TFA_output/analysis/TFA_corr/overall/(Spearman\|Pearson).pdf` |
| Pseudobulk correlations by cell type | `TFA_output/analysis/TFA_corr/pseudobulk/mean_(Spearman\|Pearson).pdf` |
| TFA and RNA feature plots | `plots/TFA/` |

In `plots/TFA/`, each subdirectory represents RNA or a TFA inference algorithm. For example, `500_100` denotes \(\lambda_C=500\) and \(\lambda_E=100\).

## PBMC10K Downstream TFA Analysis

### 4. Differential expression and activity analysis

We used Seurat's [FindAllMarkers](https://satijalab.org/seurat/reference/findallmarkers) to identify differentially expressed genes and differentially active TFs across cell types.

**Results:**

| Data | Output |
|---|---|
| TFAct-F (\(\lambda_C=500,\lambda_E=100\)) | `plots/DE/ibpg_2fct_500_100.txt` |
| RNA | `plots/DE/RNA.txt` |

Upregulated genes or more active TFs were selected using:

- Adjusted p-value `p_val_adj < 0.1`
- Average log2 fold change `avg_log2FC > 0`

**Filtered results:**

```text
plots/DE/up_regulated/ibpg_2fct_500_100/<cell_type>.txt
plots/DE/up_regulated/RNA/<cell_type>.txt
```

For downregulated genes or less active TFs, replace `up_regulated` with `down_regulated`.

### 5. Pseudotime differential analysis using tradeSeq

Initial Slingshot analysis across all PBMC cell types produced biologically implausible lineages. We therefore defined four lineage subsets and applied Slingshot separately to each.

Install dependencies for tradeSeq and run:

```bash
Rscript tradeseq.R
```

The script performs trajectory inference and differential analysis of TFA and RNA using [Slingshot](https://bioconductor.org/packages/slingshot/) and [tradeSeq](https://bioconductor.org/packages/tradeSeq/).

To compare the top 10 genes or TFs identified for each lineage, run:

```bash
Rscript tradeseq_results_jaccard.R
```

This evaluates the overlap between top-ranked findings across inference methods using Jaccard similarity.