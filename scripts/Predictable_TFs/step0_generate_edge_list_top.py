#pwd=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data
# or at the directory scMTNI_data. Make sure scMTNI_sourcedata is a subdir and that the snakemake script for generating outputs have completed

import numpy as np
import multiprocessing

top=10000 # about this number of edges will be output, as multiple edges could have the same confidence
thr = 1e-10 # for general thresholding. For thresholding kl nmf, one needs to raise the threshold value to 0.5/A.shape[1]
# TFA_output/ibpg_2fct/processed_Buenrostro_Hematopoiesis/cluster1/rep_0/A.txt

# note: regularized_nca is already in the edge-list format
nct_ss_algo = ['HALS-p', 'ibpg']
ct_ss_algo = ['ibpg_2fct_2000_1000', 'ibpg_2fct_2000_100', 'ibpg_2fct_100_1000', 'ibpg_2fct', 'kl_nmf_2fct', 'kl_nmf_2fct_2000_1000']
stability_selection_algo = nct_ss_algo + ct_ss_algo
confidence_level_algo = ['NetREX_CF']
datasets=['processed_FBS_reprogramming', 'processed_Buenrostro_Hematopoiesis', 'processed_Ranzoni_fetal_hematopoiesis_fine_tree']
FBS_clusters = ["cluster1","cluster2","cluster3","cluster4","cluster5","cluster7","cluster8"]

# Buenrostro_Hematopoiesis 
BH_clusters = ["cluster1","cluster2","cluster3","cluster6","cluster7","cluster8", "cluster9", "cluster10"]

# Ranzoni_fetal_hematopoiesis/coarse_tree
RFH_fine_clusters = ["B-cells","Erythroid-cells","GPs","Granulocytes","HSC-MPPs","LMPs","Mast-cells","Megakaryocytes","MEMPs","Monocytes","pDCs"]

ds_cluster_map = {"processed_FBS_reprogramming": FBS_clusters, "processed_Buenrostro_Hematopoiesis": BH_clusters, "processed_Ranzoni_fetal_hematopoiesis_fine_tree": RFH_fine_clusters}

mouse_gold_standards = ['mESC_chipunion_KDUnion_intersect', 'mESC_KDUnion', 'mESC_chipunion']
mouse_gold_standards_paths = ['scMTNI_sourcedata/Goldstandards/mouse_goldstandard/' + _ + '.txt' for _ in mouse_gold_standards]

# human. Work for both ds.
BH_gold_standards = ['Cus_KO_0.01_Cus_ChIP_2_intersect.gs_network', 'Cus_ChIP_2.gs_network', 'Cus_KO_0.01.gs_network', 'Cus_KO_0.01_union_intersect.gs_network', 'UniBind.gs_network']
BH_gold_standards_paths = ['scMTNI_sourcedata/Goldstandards/human_goldstandard/' + _ + '.txt' for _ in BH_gold_standards]

gs_dict={"processed_FBS_reprogramming": mouse_gold_standards, "processed_Buenrostro_Hematopoiesis": BH_gold_standards, "processed_Ranzoni_fetal_hematopoiesis_fine_tree": BH_gold_standards}
mouse_path_dict={mouse_gold_standards[i]:mouse_gold_standards_paths[i] for i in range(len(mouse_gold_standards))}
BH_path_dict = {BH_gold_standards[i]:BH_gold_standards_paths[i] for i in range(len(BH_gold_standards))}
gs_path_dict={**mouse_path_dict, **BH_path_dict}

nca_rep=20

def load_and_thr(cur_algo, ds, c, r):
    tmp = np.loadtxt('TFA_output/'+cur_algo+'/'+ds+'/'+c+'/rep_'+str(r)+'/A.txt')
    if 'kl_nmf' in cur_algo:
        start_idx=np.count_nonzero(tmp)-min(top, np.count_nonzero(tmp)/2) # just to make enough edges
        tmp_top_thr=np.partition(tmp[np.nonzero(tmp)], start_idx)[start_idx]
        tmp_top = (tmp>=tmp_top_thr) * tmp 
        return (np.abs(tmp_top) > 0).astype(int)
    return (np.abs(tmp) > thr).astype(int)

def get_confidence_edges(cur_algo, ds, c, num_rep):
    with multiprocessing.Pool() as pool_obj:
        A_mat = pool_obj.starmap(load_and_thr, [(cur_algo, ds, c, r) for r in range(num_rep)])
        A_mat = sum(A_mat, start=np.zeros_like(A_mat[0])) / num_rep
    return A_mat

# all other algos pick top edges by absolute magnitude
for ds in datasets:
    for c in ds_cluster_map[ds]:
        for cur_algo in stability_selection_algo:
            if cur_algo in ct_ss_algo: # allow union
                genes = np.loadtxt('scMTNI_sourcedata/'+ds+'/inputdata/all_exps/all_genes.txt', dtype=str)
                tfs = np.loadtxt('scMTNI_sourcedata/'+ds+'/inputdata/allregulators.txt', dtype=str)
            else:
                genes = np.loadtxt('scMTNI_sourcedata/'+ds+'/inputdata/'+ c +'_allGenes.txt', dtype=str)
                tfs = np.loadtxt('scMTNI_sourcedata/'+ds+'/inputdata/'+ c +'_allregulators.txt', dtype=str)
            A_mat_freq = get_confidence_edges(cur_algo, ds, c, nca_rep)
            if cur_algo=="HALS-p":
                A_mat_freq = A_mat_freq.T
            start_idx = max(np.count_nonzero(A_mat_freq)-top, 0)
            A_top_thr=np.partition(A_mat_freq[np.nonzero(A_mat_freq)], start_idx)[start_idx]
            A_top = (A_mat_freq>=A_top_thr) * A_mat_freq 
            sorted_indices = sorted(np.array(np.nonzero(A_top)).T, key=lambda x: A_top[x[0], x[1]], reverse=True)
            outfile = 'TFA_output/'+cur_algo+'/'+ds+'/'+c+'/A_top'+str(top)+'.txt'
            with open(outfile, 'w') as f:
                for t,r in sorted_indices:
                    _ = f.write(tfs[r]+"\t"+genes[t]+"\t"+str(A_top[t,r])+"\n")
            for gs in gs_dict[ds]:
                gs_path=gs_path_dict[gs]
                out_compare_f='TFA_output/'+cur_algo+'/'+ds+'/'+c+'/A_top'+str(top)+'.'+gs+'.txt'
                with open(out_compare_f, 'w') as f:
                    _ = f.write("gs\t"+gs_path+"\t1\n")
                    _ = f.write("inferred\t"+outfile+"\t2\n")
            print("complete "+cur_algo+"_"+ds+"_"+c)

for ds in datasets:
    for c in ds_cluster_map[ds]:
        for cur_algo in confidence_level_algo:
            genes = np.loadtxt('scMTNI_sourcedata/'+ds+'/inputdata/'+ c +'_allGenes.txt', dtype=str)
            tfs = np.loadtxt('scMTNI_sourcedata/'+ds+'/inputdata/'+ c +'_allregulators.txt', dtype=str)
            confidence= np.loadtxt('TFA_output/'+cur_algo+'/'+ds+'/'+c+'/confidence.txt')
            A = np.loadtxt('TFA_output/'+cur_algo+'/'+ds+'/'+c+'/A.txt')
            confident_and_nonzero = confidence * (A != 0).astype(int)
            start_idx = max(np.count_nonzero(confident_and_nonzero)-top, 0)
            A_top_thr=np.partition(confident_and_nonzero[np.nonzero(confident_and_nonzero)], start_idx)[start_idx]
            A_top = (confident_and_nonzero>=A_top_thr) * confident_and_nonzero 
            sorted_indices = sorted(np.array(np.nonzero(A_top)).T, key=lambda x: A_top[x[0], x[1]], reverse=True)
            outfile = 'TFA_output/'+cur_algo+'/'+ds+'/'+c+'/A_top'+str(top)+'.txt'
            with open(outfile, 'w') as f:
                for t,r in sorted_indices:
                    _ = f.write(tfs[r]+"\t"+genes[t]+"\t"+str(A_top[t,r])+"\n")
            for gs in gs_dict[ds]:
                gs_path=gs_path_dict[gs]
                out_compare_f='TFA_output/'+cur_algo+'/'+ds+'/'+c+'/A_top'+str(top)+'.'+gs+'.txt'
                with open(out_compare_f, 'w') as f:
                    _ = f.write("gs\t"+gs_path+"\t1\n")
                    _ = f.write("inferred\t"+outfile+"\t2\n")
            print("complete "+cur_algo+"_"+ds+"_"+c)
