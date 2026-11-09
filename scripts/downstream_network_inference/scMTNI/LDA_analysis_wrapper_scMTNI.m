
%% supply in the calling space
%cells=importdata('/mnt/dv/wid/projects7/Roy-singlecell2/suvo_work/scMTNI_brandlab_human_Combined/inputs_for_GRN_without_TFA/celltype_order.txt')
%indir='/mnt/dv/wid/projects7/Roy-singlecell2/suvo_work/scMTNI_brandlab_human_Combined/inputs_for_GRN_without_TFA/results_b_4_q_1/analysis/lda_TFcellbygene/'
%cf='_cf0.7';
%k_lst=[15 22 30]


%%cf='_top2000';
prefix='_filteredlowexpression'
dat=''
binry=0
%% No. of rows in the "lda_k*__tf-topicPercell_filteredlowexpression*.pdf" plot
nrp=3
%% No. of columns in the "lda_k*__tf-topicPercell_filteredlowexpression*.pdf" plot
ncp=4

%% k = Number of topics. Please see MATLAB function fitlda().
%for k=[30,25,20,15]
for k=k_lst
	fprintf('\n\n\n starting k=%d \n\n\n',k)
    [mdl10]=LDA_analysis_v2(indir,cells,k,prefix,cf,dat,binry,nrp,ncp)
	fprintf('\n\n\n done with k=%d \n\n\n',k)
end














