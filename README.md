# TFAct: TFA with cell type lineage

TFAct is a matrix factorization-based method for inferring cell type-specific gene regulatory interactions and estimating the latent activity of transcription factors.

## Installation

### Option 1:  (without conda) 
**step 1**: install [gsl](https://www.gnu.org/software/gsl/)

```bash
# can install a version >= 2.6
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.7.tar.gz
tar -zxvf gsl-2.7.tar.gz
cd gsl-2.7

# a directory where you would like gsl to be installed
mkdir /home/yourname/gsl
./configure --prefix=/home/yourname/gsl

# install
make
make check
# if passed, proceed:
make install
```
**step 2**: cd into the TFAct directory (this repo) and replace the definition of CONDA_PREFIX in **Makefile** with the prefix /home/yourname/gsl specified earlier:
```bash
CONDA_PREFIX = /home/yourname/gsl
```

**step 3**: make
```bash
make
```
### Option 2:  (with conda) 
**step 1**: [create](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) a conda environment with name or prefix. With the environment activated, do
```bash
conda install -c conda-forge gsl
```
**step 2**: replace the definition of CONDA_PREFIX in **Makefile** with the conda prefix. If the environment was created by name, the conda prefix can be found by the environment variable **$CONDA_PREFIX**.

**step 3**: make
```bash
make
```

## Usage
Before running, ensure that you have specified the LD_LIBRARY_PATH 
```bash
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CONDA_PREFIX}/lib
```
To apply TFAct to the example files with 10 random inits:
```bash
./tetraNMF \
  -d example/in/all_exps_differentiation_0.2_k20/E_all_exps.txt \
  -r example/in/all_exps_differentiation_0.2_k20/tflist.txt \
  -g example/in/all_exps_differentiation_0.2_k20/all_genes.txt \
  -p example/in/all_exps_differentiation_0.2_k20/all_priors.txt \
  -o example/out/all_exps_differentiation_0.2_k20 \
  --max_iter 300 --init_m random \
  --num_rep 10 -t ibpg_2fct --lambda_p 0.1 --lambda_Cminus 0.1 \
  --tree example/in/all_exps_differentiation_0.2_k20/cell_lineage.txt \
  -l 500 --lambda_E 100
```

## Output

rep_[0-9]+/${cell_type}/A.txt: GRN matrices for ${cell_type} for the random init.

rep_[0-9]+/${cell_type}/P.txt: TFA matrices for ${cell_type} for the random init.
