## Downstream GRN Inference Using TFA

We evaluated whether TFA estimates from different algorithms improve gene regulatory network (GRN) inference using **Inferelator** and **MERLIN-P** on the scMTNI datasets. TFA estimates were obtained from the preceding experiments.

### 1. Install Inferelator

```bash
mamba create --prefix ./env/inferelator3 python=3.10 -y
mamba activate ./env/inferelator3
python -m pip install inferelator
```

### 2. Configure HTCondor

We used the [Snakemake HTCondor profile](https://github.com/Snakemake-Profiles/htcondor) to submit inference jobs.

**Requirements:**

- Condor-submitted rules must specify disk, memory, and thread requirements. Rules without these resources should run locally.
- Remote workers must have access to the shared filesystem. We therefore restricted execution to WID cluster machines.
- MERLIN-P was submitted separately using `submit_files/submit_merlinp_p-5_r4.sub`.

Workflow configuration and input preparation are described in `Snakefile` and `rules/`.

### 3. Run MERLIN-P

**Prerequisite:** Generate MERLIN-P input files as described in the comments at the beginning of `Snakefile`.

Run from the project directory on a Condor submit node:

```bash
bash make_queue_args.sh
condor_submit submit_files/submit_merlinp_p-5_r4.sub
```

**Important:** Run one block at a time in `make_queue_args.sh`, since each block corresponds to a dataset. Block 4 contains the expression-only configuration and can be run first.

#### Organize outputs and rerun failed jobs

After Condor jobs complete:

```bash
bash ../../scripts/downstream_network_inference/organize_merlin_condor_output.sh
```

This script identifies invalid networks (including rows with a number of fields other than three) and jobs that produced no valid output.

Failed-job logs are saved under:

```text
GRN_outputs/failed_merlin_jobs/logs/
GRN_outputs/failed_merlin_jobs/logs_not_outputting/
```

Resubmission lists are generated under:

```text
GRN_outputs/failed_merlin_jobs/*queue_args.txt
```

These include dataset-specific files for:

- `processed_FBS_reprogramming`
- `processed_Buenrostro_Hematopoiesis`
- `processed_Ranzoni_fetal_hematopoiesis_fine_tree`

To rerun failed jobs, copy the relevant non-empty queue argument files into `queue_args.txt` and resubmit:

```bash
cp GRN_outputs/failed_merlin_jobs/<failed_queue_args_file> queue_args.txt
condor_submit submit_files/submit_merlinp_p-5_r4.sub
```

Repeat output organization and resubmission until all required networks are generated.

**Troubleshooting:** Some jobs terminated with exit code 134 (`SIGABRT`). In `run_merlinp_condor.sh`, cleanup handling can be extended from:

```bash
trap clean_up SIGTERM
```

to:

```bash
trap clean_up SIGTERM SIGABRT
```

Some observed failures were potentially related to cluster server issues rather than the inference algorithm.

### 4. Run Inferelator

**Prerequisites:** Complete input preparation according to `Snakefile` and configure the HTCondor profile. Ensure MERLIN-P targets are commented out in the `all` rule when running Inferelator alone.

```bash
snakemake --profile htcondor --use-conda -k --rerun-incomplete
```

The workflow is defined primarily in `rules/inferelator.smk`.

We evaluated three types of regulator features:

1. TFA values inferred by different TFA algorithms.
2. Inferelator's internally estimated TFA.
3. `no_TFA`, using regulator expression instead of inferred activity.

Inferelator is applied to cell-type-specific subsamples independently, using the same regulator sets as the TFA inference algorithms.

For grouping Condor jobs by Snakemake rule, the profile's `grid-submit.py` was modified (see the local HTCondor profile configuration).

### 5. Inferelator evaluation

Evaluation includes AUPR, regulator-specific AUPR, F-scores, and predictable TF counts.

**Output directories** (relative to `GRN_outputs/algos/analysis/`):

| Metric | Directory |
|---|---|
| Overall AUPR | `inferelator_AUPR/` |
| Regulator-specific mean/median AUPR | `inferelator_AUPR_TF/` |
| Predictable TFs | `inferelator_predictable_TFs/` |
| F-scores | `inferelator_Fscore/` |

AUPR and F-score results include heatmaps and critical difference plots. We followed [Demšar (2006)](https://www.jmlr.org/papers/volume7/demsar06a/demsar06a.pdf) and used the corresponding [R implementation](https://journal.r-project.org/archive/2016/RJ-2016-017/index.html) for comparisons across datasets.

Generate plots:

```bash
Rscript ../../scripts/downstream_network_inference/plot_AUPR_inferelator_scMTNIdata.R

Rscript ../../scripts/downstream_network_inference/plot_AUPR_TF_inferelator_scMTNIdata.R
```

**Representative outputs:**

```text
inferelator_AUPR/
    AUPR_<gold_standard>_<dataset>_critical_difference.pdf
    <dataset>_<gold_standard>.pdf

inferelator_Fscore/
    Fscore_top_<500|1000|2000>_<gold_standard>_<dataset>_critical_difference.pdf
    <dataset>_<gold_standard>_<500|1000|2000>.pdf
```

#### Predictable TFs

Requires Python with the necessary analysis libraries and GSL shared libraries.

```bash
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/envs/IBPG/lib"

top_nums=(5000 10000 30000 50000 100000)

for top in "${top_nums[@]}"; do
    bash ../../scripts/downstream_network_inference/inferelator_predictable_TFs/step1_compute_predictable_TFs.sh "$top" &
done
wait

for top in "${top_nums[@]}"; do
    python ../../scripts/downstream_network_inference/inferelator_predictable_TFs/step2_heatmap_plot_predictable_TFs.py -t "$top" &
done
wait
```

### 6. Compare Inferelator GRNs with TFA-derived A matrices

We compared the mean and median TF-AUPR of Inferelator-inferred networks against the TF-target connectivity matrices (A matrices) learned during TFA inference.

First, generate TF-AUPR for the TFA inference algorithms from the `scMTNI_data` experiment directory:

```bash
Rscript ../../scripts/plot_AUPR_TF_scMTNI.R
```

Also complete the Inferelator TF-AUPR analysis in Step 5.

Then, from the downstream GRN inference directory:

```bash
Rscript ../../scripts/downstream_network_inference/AUPR_TF_btw_A_matrix_and_inferelator_GRN.R
```

The script performs one-sided Wilcoxon signed-rank tests, with the alternative hypothesis that Inferelator improves regulator-specific AUPR.

**Output directory:**

```text
GRN_outputs/algos/analysis/inferelator_AUPR_TF/comparison_w_Amat/
```

**Outputs:**

```text
(mean|median)_<dataset>_<gold_standard>_wilcox_signed_rank.txt
(mean|median)_<gold_standard>_wilcox_signed_rank.txt
(mean|median)_<dataset>_wilcox_signed_rank.txt
```

These summarize comparisons by dataset and gold standard, across gold standards, and across the two human datasets, respectively (see the R script for the exact naming conventions).

We observed that methods other than TFAct-KL and SupirFactor generally did not exhibit an improvement in regulator-specific AUPR after Inferelator inference.

### 7. MERLIN-P evaluation

MERLIN-P can use regulators without estimated TFA. We therefore used the scMTNI regulator sets for each cell type and compared two ways of incorporating TFA:

- **Replaced:** Replace TF expression with estimated TFA when available.
- **Concatenated:** Append TFA as additional features alongside gene expression.

For regulators with estimated TFA, their expression features were retained when constructing the regulator set. Regulators not expressed in at least 50 cells across any cell cluster were excluded.

For concatenated inputs, we extract three networks from each inferred consensus network:

1. `gene_only`: mRNA -> mRNA edges.
2. `TFA_only`: TFA -> mRNA edges.
3. `max_conf`: For each regulator-target pair, retain the larger confidence score across gene-only and TFA-only edges.

#### Process MERLIN-P outputs

After completing and organizing Condor jobs, uncomment `rules.merlin_output.input` in the `all` rule of `Snakefile`.

```bash
snakemake --cores 30 -k --rerun-incomplete --use-conda
```

This generates the processed networks and evaluation inputs.

#### Generate evaluation plots

```bash
Rscript ../../scripts/downstream_network_inference/plot_AUPR_merlin_scMTNIdata.R

Rscript ../../scripts/downstream_network_inference/plot_AUPR_TF_merlin_scMTNIdata.R
```

**Output directories** (under `GRN_outputs/algos/analysis/`):

```text
merlin_replaced_AUPR/
merlin_replaced_AUPR_TF/
merlin_replaced_Fscore/
merlin_replaced_predictable_TFs/

merlin_concat_AUPR/
merlin_concat_AUPR_TF/
merlin_concat_Fscore/
merlin_concat_predictable_TFs/
```

These contain the corresponding evaluation summaries, plots, and method rankings.

### 8. MERLIN-P predictable TFs

Set the GSL library path if needed:

```bash
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/envs/IBPG/lib"
```

Compute predictable TFs for different edge cutoffs:

```bash
top_nums=(5000 10000 30000 50000 100000)

for top in "${top_nums[@]}"; do
    bash ../../scripts/downstream_network_inference/merlin_predictable_TFs/step1_compute_predictable_TFs.sh "$top" &
done
wait
```

#### Replaced TFA

```bash
for top in "${top_nums[@]}"; do
    python ../../scripts/downstream_network_inference/merlin_predictable_TFs/step2_heatmap_plot_predictable_TFs.py \
        -t "$top" \
        -i "merlin_predictable_tf/c5_p-5_r4/merlin_predictable_tf_nums_top${top}.txt" \
        -o "GRN_outputs/algos/analysis/merlin_replaced_predictable_TFs/top_${top}" &
done
wait
```

#### Concatenated TFA

Generate predictable TF heatmaps separately for gene-only, TFA-only, and maximum-confidence networks:

```bash
for etype in gene_only TFA_only max_conf; do
    for top in "${top_nums[@]}"; do
        python ../../scripts/downstream_network_inference/merlin_predictable_TFs/step2_heatmap_plot_predictable_TFs.py \
            -t "$top" \
            -i "merlin_predictable_tf/c5_p-5_r4/concat/${etype}_merlin_predictable_tf_nums_top${top}.txt" \
            -o "GRN_outputs/algos/analysis/merlin_concat_predictable_TFs/${etype}_top_${top}" &
    done
done
wait
```

All commands assume execution from the `GRN_inference_w_TFA` experiment directory, with analysis scripts available under `../../scripts/downstream_network_inference/`.