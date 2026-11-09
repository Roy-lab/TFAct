function [A,S] = ncanec_quad_sub_nocvx_colwise( X, A_top, method, tls_scaling )
% [A,S] = ncanec_quad_sub_nocvx_colwise( X, A_top ) computes the recovered 
% connectivity matrix A and the recovered source signals S given 
% observations X.  The network size is inferred from the input arguments.


if ( nargin < 3 )
    method = 'tls';
end

if ( nargin < 4 )
    tls_scaling = 1e4;
end


% find N, M, K; check consistency
[N,M,K] = helper_get_network_size(X,A_top);

% find zero locations in A from A_top
c = 0.001; % "zero" threshold
I = find( A_top  < c );
J = find( A_top >= c );


% run NCA-NEC algorithm

[U,E,V] = svd(X);

U_S = U(:,1:M); % signal subspace
U_0 = U(:,M+1:end); % noise subspace

A = zeros(N,M);

for j = 1:M % for each col of A, called j
    a_j = A_top(:,j);
    
    H = U_0*U_0';
    
    % a_j(I) == 0
    I_j = find( a_j == 0 );
    aeq0 = zeros(length(I_j),N);
    beq0 = zeros(length(I_j),1);
    for i_j_index = 1:length(I_j)
        aeq0(i_j_index,I_j(i_j_index)) = 1;
    end
    
    % sum( a_j ) = || a_j ||_0 (from A_topology)
    aeq1 = ones(1,N);
    beq1 = sum( A_top(:,j) > 0 );
    
    % build constraints
    aeq = [aeq0; aeq1];
    beq = [beq0; beq1];
    
    a_j_opt = quadprog(H,[],[],[],aeq,beq,zeros(N,1),[],a_j);
    
    A(:,j) = a_j_opt;
    
end

% Reconstruct 'source' signals S
S = reconstruct_sources(A,X,method,tls_scaling);
