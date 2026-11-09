import argparse
import numpy as np
from inferelator import inferelator_workflow, MPControl

def run_inferelator(args):
    worker = inferelator_workflow(regression="bbsr", workflow="tfa")
# worker.set_file_paths(output_dir="./output_inferelator3",
#                   expression_matrix_file="/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data/scMTNI_sourcedata/processed_FBS_reprogramming/inputdata/subsamples/cluster1/dataset.0.txt",
#                   tf_names_file="/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data/scMTNI_sourcedata/processed_FBS_reprogramming/inputdata/cluster1_allregulators.txt",
#                   meta_data_file=None,
#                   priors_file="/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data/scMTNI_sourcedata/processed_FBS_reprogramming/inputdata/inferelator.cluster1_network.tsv.gz",
#                   gold_standard_file="/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data/scMTNI_sourcedata/processed_FBS_reprogramming/inputdata/inferelator.cluster1_network.tsv.gz")
    worker.set_file_paths(output_dir=args.outdir,
                      expression_matrix_file=args.exp,
                      tf_names_file=args.tf_lst,
                      meta_data_file=None,
                      priors_file=args.prior_A,
                      gold_standard_file=args.gs_A)
    worker.set_file_properties(expression_matrix_columns_are_genes=False)
    worker.set_run_parameters(num_bootstraps=5, random_seed=42)
    if args.tfa is not None:
        # worker.set_tfa(tfa_driver=True, tfa_input_file='test/test_TFA.tsv.gz', tfa_input_file_type='tsv')
        worker.set_tfa(tfa_driver=True, tfa_input_file=args.tfa, tfa_input_file_type='tsv')
    if args.no_TFA:
        worker.set_tfa(tfa_driver=False)
    network_result = worker.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp", help="a tsv expression matrix.", type=str, required=True)
    parser.add_argument("--tf_lst", help="a column file where each line is the name of a regulator.", type=str, required=True)
    parser.add_argument("--prior_A", help="The prior GRN.", type=str, required=True)
    parser.add_argument("--gs_A", help="The gold standard GRN. We don't want inferelator to have access to gold standard during training, so can just supply the prior.", type=str, required=True)
    parser.add_argument("--outdir", help="output directory.", type=str, required=True)
    parser.add_argument("--tfa", help="sample by tfa tsv file path", type=str)
    parser.add_argument("--no_TFA", help='Use TF expression instead', action='store_true')
    parser.add_argument("--threads", help='Number of threads to use', type=int, required=True)

    args = parser.parse_args()
    MPControl.set_multiprocess_engine("multiprocessing")
    MPControl.client.processes = args.threads
    MPControl.connect()
    run_inferelator(args)