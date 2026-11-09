function run_ncanec_quad_sub_nocvx_colwise(X_path, prior_path, tf_path, A_out, P_out)

X_in = importdata(X_path);
A_in = importdata(prior_path);
tf=importdata(tf_path);


gene=X_in.textdata;

num_rows = length(gene);
num_cols = length(tf);
A = zeros(num_rows, num_cols);

gene_csr = A_in.textdata(:, 2);  
tf_csr = A_in.textdata(:, 1);  

values = A_in.data;

[~, row_indices] = ismember(gene_csr, gene);
[~, col_indices] = ismember(tf_csr, tf);

for i = 1:length(values)
    row_idx = row_indices(i);
    col_idx = col_indices(i);
    A(row_idx, col_idx) = values(i);
end


X = X_in.data;

[Ae,Se] = ncanec_quad_sub_nocvx_colwise(X, A);
dlmwrite(A_out, Ae, 'delimiter', ' ');
dlmwrite(P_out, Se, 'delimiter', ' ');
