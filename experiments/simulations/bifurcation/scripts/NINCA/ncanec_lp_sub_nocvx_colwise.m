function [A,S] = ncanec_lp_sub_nocvx_colwise( X, A_top, method, tls_scaling )
% [A,S] = ncanec_lp_sub_nocvx_colwise( X, A_top ) computes the recovered 
% connectivity matrix A and the recovered source signals S given 
% observations X.  The network size is inferred from the input arguments.


if ( nargin < 3 )
    method = 'inverse';
end

if ( nargin < 4 )
    tls_scaling = 1;
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

% build up params for linprog()
%     X=LINPROG(f,A,b) attempts to solve the linear programming problem:
%          
%              min f'*x    subject to:   A*x <= b 
%               x
%  
%     X=LINPROG(f,A,b,Aeq,beq) solves the problem above while additionally
%     satisfying the equality constraints Aeq*x = beq.
%     
%     X=LINPROG(f,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
%     bounds on the design variables, X, so that the solution is in
%     the range LB <= X <= UB.  Use empty matrices for LB and UB
%     if no bounds exist. Set LB(i) = -Inf if X(i) is unbounded below; 
%     set UB(i) = Inf if X(i) is unbounded above.
%  
%     X=LINPROG(f,A,b,Aeq,beq,LB,UB,X0) sets the starting point to X0.  This
%     option is only available with the active-set algorithm.  The default
%     interior point algorithm will ignore any non-empty starting point.

% Solve column-wise.  
% Note: optimization variable will be x = [t; a_j]


A = zeros(N,M);

f = [1; zeros(N,1)];
bc = zeros(2*(N-M),1);

for j = 1:M % for each col of A, called j
    a_j = A_top(:,j);
    
    % -t < U_0(:,i)'*a_j < t constraints, i = 1:(N-M)
    Ac = zeros(2*(N-M), 1+N); % rows in Ac = # of ineq. constraints
    
    % U_0(:,i)'*a_j - t < 0
    c_index = 1;
    for i = 1:N-M % row of U_0'
        Ac(i,:) = [-1 U_0(:,i)'];
    end
    
    % - U_0(:,i)'*a_j - t < 0
    for i = 1:N-M
        Ac((N-M)+i,:) = [-1 -U_0(:,i)'];
    end
    
    % a_j(I) == 0
    I_j = find( a_j == 0 );
    aeq0 = zeros(length(I_j),1+N);
    beq0 = zeros(length(I_j),1);
    for i_j_index = 1:length(I_j)
        aeq0(i_j_index,1+I_j(i_j_index)) = 1;
    end
    
    % sum( a_j ) = || a_j ||_0 (from A_topology)
    aeq1 = [0 ones(1,N)];
    beq1 = sum( A_top(:,j) > 0 );
    
    % build constraints
    aeq = [aeq0; aeq1];
    beq = [beq0; beq1];
    
    a_j_opt = linprog(f,Ac,bc,aeq,beq,zeros(1+N,1),[],a_j); % need 2022b, or omit a_j. Any way a_j will not be used
    
    A(:,j) = a_j_opt(2:end);
    
end

% Reconstruct 'source' signals S
S = reconstruct_sources(A,X,method,tls_scaling);

end

