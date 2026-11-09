function wrapperNMF_simple(in_file, k_range, num_rand_init, title, out_dir)
%% this code is adapted from /mnt/dv/wid/projects7/Roy-plants/Sonar_project/saptarshi_work/scripts/wrapperNMF_simple_MtRafael.m
%% Usage:
%% module load matlab-r2017b
%% matlab -nodisplay -nodesktop -r "run wrapperNMF_simple.m"
%%
%% Goal:
%% Apply simple NMF on a given data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin: Set input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data file
% in_file = '/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/simulation/072024_cell_lineage_formulation_w_lda_noise/inputtoTFA/0.2/k10/Anc2_prior.txt';

%% Range of k values
% k_range = [10];


%% Number of random initializations to be performed for each value of k.
%% It is necessary for checking the stability of NMF.
% num_rand_init = 10;

% title = 'Anc2_k10_prior';

% out_dir = 'test3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End: Set input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import (tfs x genes) edge file, where each line is an edge
% d = importdata(in_file,'\t');
% data = d.data(:,:);
tdata = readtable(in_file, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);

%% Retrieve tf and gene names.
tf_names = unique(tdata.Var1);
gene_names = unique(tdata.Var2);

m = length(tf_names);
n = length(gene_names);
parentMap = containers.Map(tf_names, 1:m);
childMap = containers.Map(gene_names, 1:n);

d = zeros(m, n);

for i = 1:size(tdata, 1)
        parent = string(tdata.Var1(i));
        child = string(tdata.Var2(i));
        r = parentMap(parent);
        c = childMap(child);
        d(r, c) = 1;
end

adj_data=d;

for rand_init_no = 1:num_rand_init
   for k_idx = 1:length(k_range)

        k_val = k_range(k_idx);

	%% Create the output sub-dir for the current random initialization and the current value of k
        out_subdir = sprintf('%s/K%d/I%02d', out_dir, k_val, rand_init_no);
        mkdir(out_subdir);

	%% Run NMF not in parallel
        [U, V] = nnmf(adj_data, k_val, 'Options', statset('UseParallel',false));

	%% Save the NMF factors inside a file
        outs_depthnormalize{1} = U;
        outs_depthnormalize{2} = V;
        save(sprintf('%s/outs_depthnormalize', out_subdir), 'outs_depthnormalize');

	%% Calculate barcode and feature clusters separately.
	%% Save the cluster assignments in two seprate files.
	%% Also, generate a heatmap of NMF factorization for each random initialization in each value of k.
        addpath('/mnt/dv/wid/projects5/Roy-singlecell/sr_work/nmf_analysis');   % makefigure_newvis.m
        % [rid, cid, rorder, corder] = makefigure_newvis(full(u),full(v),alldata,tname,nparts,0,3,goodregionnames,tf_names,outsubdir);
	[rid, cid, rorder, corder] = makefigure_newvis(U, V, adj_data, title, k_val, 0, 3, gene_names, tf_names, out_subdir);
   end
end

