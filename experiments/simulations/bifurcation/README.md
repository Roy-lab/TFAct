# Simulated bifurcation instructions


This pipeline generates ground-truth regulatory networks, noisy prior networks, and simulated transcription factor activity (TFA) and gene expression data.

## Step 1: Generate ground-truth networks

We simulate two edge gaining probabilities: **0.2 (denser)** and **0.1 (Sparser)**, corresponding to the bifurcation scheme in the manuscript. The edge retaining probability is 0.8.

```bash
indir=../speciesspecsim_sz
mkdir -p testoutput

./${indir}/simulateNw ${indir}/testinput/net300.tsv species_prob_pgain0.2.txt ${indir}/testinput/tflist.txt testoutput/0.2/ 0.2

./${indir}/simulateNw ${indir}/testinput/net300.tsv species_prob_pgain0.1.txt ${indir}/testinput/tflist.txt testoutput/0.1/ 0.2
```

The final parameter (`0.2`) approximately controls the fraction of edges from the complete TF-target graph (net300.tsv) retained in the root network. Larger values produce denser networks. The generated networks are treated as **ground truth**.

## Step 2: Generate noisy prior networks using LDA

LDA is used to introduce noise into the ground-truth networks to make noisy prior TF-target connectivity.

- **K = 10:** Higher noise
- **K = 20:** Lower noise
- **Sparsity parameter = 1.5:** Larger values produce denser networks.

```bash
module load matlab-r2023a # require matlab

simdir=$(pwd)

for density in 0.2 0.1; do
    mkdir -p inputtoTFA/${density}/k10
    bash get_simulated_prior_networks.sh \
        ${simdir}/testoutput/${density} \
        inputtoTFA/${density}/k10 \
        10 1.5
    
    mkdir -p inputtoTFA/${density}/k20
    bash get_simulated_prior_networks.sh \
        ${simdir}/testoutput/${density} \
        inputtoTFA/${density}/k20 \
        20 1.6
done
```

**Output:** `inputtoTFA/<density>/k<K>/*_prior.txt`

Edges sampled once receive confidence `0.8`; edges sampled more than once receive confidence `1`.

## Step 3: Simulate TFA and gene expression

Using the ground-truth networks from Step 1, we simulate TFA and gene expression.

TFA is sampled using a truncated normal distribution:

```python
truncnorm(-0.8, 1, loc=1, scale=1)
```

This upperbounds TFA by `2`.

**Dependencies:** Python with NumPy and SciPy.

```bash
simdir=$(pwd)

for density in 0.2 0.1; do
    for k in 10 20; do
        bash get_simulated_TFA_and_expression.sh \
            ${simdir}/testoutput/${density} \
            inputtoTFA/${density}/k${k}
    done
done
```

### Output structure

```text
inputtoTFA/
└── <density>/             # 0.2 or 0.4
    └── k<K>/              # k10 or k20
        └── <celltype>/    # Calb, Sbay, Anc2, Cgla, right, Spom, Anc5
            └── rep<N>/
                ├── E.txt
                ├── E_0.6.txt
                ├── P.txt
                ├── cell_celltype_assignments.txt
                └── tf_celltype_assignments.txt
```

| File | Description |
|---|---|
| `E.txt` | Simulated gene expression matrix |
| `E_0.X.txt` | Expression matrix with dropout probability X (e.g., 0.9 means 90% dropout for positive values) |
| `P.txt` | Ground-truth TFA matrix |
| `cell_celltype_assignments.txt` | Ground-truth cell type assignment for each cell |
| `tf_celltype_assignments.txt` | Cell type in which each TF is designated active |

The original TF list is located at `${indir}/testinput/tflist.txt` from Step 1.


## Step 4: Run TFA Inference Algorithms

### 4.1 Generate inputs for TFAct

Prepare the cell lineage files and inputs for TFAct (IBPG-2FCT) across network densities and LDA noise levels.

```bash
simdir=$(pwd)

for density in 0.1 0.2; do
    for k in 10 20; do
        outdir=inputtoTFA/${density}/k${k}/all_exps
        mkdir -p "${outdir}"

        cut -f1,2 \
            species_prob_pgain${density}.txt |
            tail -n +2 |
            awk -F'\t' '{printf "%s\t%s\t1\n", $1, $2}' \
            > "${outdir}/cell_lineage.txt"

        # it is the script dir within the directory
        bash scripts/IBPG_2FCT/make_inputs.sh \
            "${simdir}/testoutput/${density}" \
            "${simdir}/inputtoTFA/${density}/k${k}" \
            "${simdir}/${outdir}"
    done
done
```

**Output:** `inputtoTFA/<density>/k<K>/all_exps/`

The `cell_lineage.txt` files describe the lineage relationships required by IBPG-2FCT.

**Note:** This experiment uses network densities `0.1` and `0.2` (edge addition probability), rather than `0.2` and `0.4` in simulated differentiation (which stands for edge retaining probability). Ensure the corresponding simulation inputs have been generated.

### 4.2 Run all TFA inference algorithms

All inference algorithms are executed using Snakemake.

**Dependencies:** Snakemake and Conda. Ensure that LD_LIBRARY_PATH is set to have gsl library. For environments required by individual algorithms, check each `rules/<algorithm>.smk`. Some require matlab (module load ...). For others that require conda environment, install the environments and specify the paths under the conda directive for each rule.

Before running snakemake, modify the config.yml file to specify genedir and datadir following the instructions, and specify tflist, outdir, and scriptdir to be absolute paths. Complete dir_path for pmf_grn.

```bash
snakemake -k --cores 70 --use-conda --rerun-incomplete
```

The command executes the workflow with up to 70 CPU cores, uses Conda environments, continues independent jobs after failures (`-k`), and reruns incomplete jobs (`--rerun-incomplete`).

Run the command from the experiment directory containing the `Snakefile`.

## Step 5: Visualize results
We use R to visualize the benchmarking results. Install dependencies as needed.

```bash
## GRN
Rscript scripts/plot_Fscores.R # F scores
## TFA
Rscript scripts/plot_TFA_corr.R # correlations to ground-truth
```