## scMTNI Datasets

Steps for getting benchmarking results for the real datasets without downstream GRN inference.

### 1. Data preparation

We used the published datasets from the [scMTNI study](https://www.nature.com/articles/s41467-023-38637-9#Sec35).

Because TFA inference requires prior TF-target connectivity, TFs without prior edges were excluded. This reduces the number of regulators to approximately 300–500 per dataset.

**Note:** For all NMF-based algorithms, we used random initialization instead of NNDSVD, as NNDSVD produces sparse initializations that may lead to NaN values during optimization.

### 2. Run TFA inference algorithms

The Snakemake workflow requires the following files:

```text
scMTNI_data/
├── Snakefile
├── config.yml
└── rules/
    └── *.smk
```

Set `datadir` and `scriptdir` to their corresponding absolute paths in `config.yml`.

Before running, configure the GSL shared library path (replace the path below to your gsl library which contains libgsl.so.0):

```bash
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/mnt/dv/wid/projects2/Roy-common/programs/thirdparty/gsl-2.6/lib"
```

Run Snakemake from the `scMTNI_data` directory (ensure that the matlab modules, paths to conda environment in `rules/<algorithm>.smk`, and config.yml are correctly specified as in the simulation README.md):

```bash
snakemake --cores 50 --use-conda -k --rerun-incomplete
```

**Key workflow rules:**

| Rule | Description |
|---|---|
| `*_reformatting_*` | Preprocess and reformat downloaded datasets |
| `TFA_2fct_reformat_lineage_tree` | Prepare lineage tree information |
| `TFA_2fct_make_input` | Generate inputs for TFAct-F (`IBPG_2fct`) and TFA-KL (`kl_nmf_2fct`) |
| `run_${algo}_TFA_inference` | Execute TFA inference algorithms |
| `stability_selection_get_mean` | Average replicates with sufficiently good model fit based on loss |

Algorithm-specific rules are defined in `rules/${algo}.smk`.

### 3. Plot evaluation results

Run the following scripts from the `scMTNI_data` directory:

```bash
Rscript ../../scripts/plot_AUPR_TF_scMTNI.R
Rscript ../../scripts/plot_Fscores_scMTNI.R
Rscript ../../scripts/plot_AUPR_scMTNI.R
```

**Note on stability selection:**

- F-scores are computed by selecting top-ranked edges based on their absolute magnitudes, rather than their frequencies across runs.
- Predictable TFs are evaluated using edge occurrence frequencies across runs.

## Predictable TF Analysis

### 4. Compute predictable TFs

After the TFA inference algorithms finish, run the following scripts from the `scMTNI_data` directory.

Ensure `LD_LIBRARY_PATH` includes the GSL library directory as specified above.

```bash
python ../../scripts/Predictable_TFs/step0_generate_edge_list_top.py

bash ../../scripts/Predictable_TFs/step1_compute_predictable_TFs.sh

python ../../scripts/Predictable_TFs/step2_heatmap_plot_predictable_TFs.py
```

**Pipeline steps:**

1. `step0_generate_edge_list_top.py`: Generate ranked edge lists for enrichment analysis.
2. `step1_compute_predictable_TFs.sh`: Compute predictable TFs.
3. `step2_heatmap_plot_predictable_TFs.py`: Visualize predictable TF results as heatmaps.

**Important:** Currently, the AUPR computation step also generates intermediate files required for predictable TF analysis. Consequently, algorithms other than stability-selection methods and NetREX_CF (which provides edge confidence scores) must also pass through this step to generate the required inputs, even when AUPR is not directly evaluated. These two operations could be decoupled in a future implementation.