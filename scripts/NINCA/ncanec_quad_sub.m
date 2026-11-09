function [A,S] = ncanec_quad_sub( X, A_top );
% [A,S] = ncanec_quad_sub( X, A_top ) computes the recovered connectivity matrix A
% and the recovered source signals S given observations X.  The network
% size is inferred from the input arguments.

% find N, M, K; check consistency
[N,M,K] = helper_get_network_size(X,A_top);

% find zero locations in A from A_top
c = 0.01; % "zero" threshold
I = find( A_top  < c );
J = find( A_top >= c );


% run NCA-NEC algorithm

[U,E,V] = svd(X);

U_S = U(:,1:M); % signal subspace
U_0 = U(:,M+1:end); % noise subspace

cvx_begin quiet
    variable A(N,M);
    
    minimize( norm(U_0'*A,'fro') );
    
    subject to
        % zeros
        for i = I
            A(i) == 0;
        end
        
        % non-negativity
        for j = J
            A(j) >= c;
        end
        
        % avoid nontrivial solution
        for c = 1:M
            sum(A(:,c)) == 1;
        end
        
cvx_end

A = full(A); % if A is sparse, convert it to a full matrix for pinv

% Reconstruct 'source' signals S
S = pinv(A)*X; 
