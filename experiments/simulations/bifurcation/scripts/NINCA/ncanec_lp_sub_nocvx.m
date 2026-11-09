function [A,S] = ncanec_lp_sub_nocvx( X, A_top, rel_noise_power, method )
% [A,S] = ncanec_lp_sub( X, A_top ) computes the recovered connectivity matrix A
% and the recovered source signals S given observations X.  The network
% size is inferred from the input arguments.

% find N, M, K; check consistency
[N,M,K] = helper_get_network_size(X,A_top);

% find zero locations in A from A_top
c = 0.001; % "zero" threshold
I = find( A_top  < c );
J = find( A_top >= c );


if ( nargin < 3 )
    rel_noise_power = 0;
end

if ( nargin < 4 )
    method = 'inverse';
end


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

% Note: optimization variable will be x = [t; vec(A)]

% min t
f = [1; zeros(N*M,1)];

% -t < U_0(:,i)'*a_j < t constraints, i = 1:(N-M), j = 1:M
Ac = zeros(2*(N-M)*M, 1+N*M); % rows in Ac = number of ineq. constraints

% U_0(:,i)'*a_j - t < 0 
c_index = 1;
for i = 1:N-M % i determines row of U_0'
    for j = 1:M % j determines col of A
        Ac(c_index,:) = [-1 zeros(1,(j-1)*N) U_0(:,i)' zeros(1,(M-j)*N)];
        c_index = c_index + 1;
    end
end

% - U_0(:,i)'*a_j - t < 0 
for i = 1:N-M % i determines row of U_0'
    for j = 1:M % j determines col of A
        Ac(c_index,:) = [-1 zeros(1,(j-1)*N) -U_0(:,i)' zeros(1,(M-j)*N)];
        c_index = c_index + 1;
    end
end

bc = zeros(2*(N-M)*M,1);

% A(I) == 0 
Aeq0 = zeros(length(I),1+N*M);
beq0 = zeros(length(I),1);
A_0 = zeros(size(A_top));
for i_index = 1:length(I)
    A_0(I(i_index)) = 1;
    Aeq0(i_index,:) = [1 vec(A_0)'];

    A_0(I(i_index)) = 0; % reset entry to zero
end

% sum(cols(A)) == 1
Aeq1 = zeros(M,1+N*M);
beq1 = ones(M,1);

for m = 1:M
    A_0 = zeros(size(A_top));
    A_0(:,m) = 1;
    Aeq1(m,:) = [1 vec(A_0)'];
end

% build Aeq, beq
Aeq = [ Aeq0; Aeq1 ]; 
beq = [ beq0; beq1 ];

% A(J) >= c
lb = zeros(1+N*M,1);
[J_r,J_c] = ind2sub(size(A_top),J);
lb(1) = -Inf;
for j_index = 1:length(J)
    lb( 1 + J_r(j_index) + N*(J_c(j_index)-1) ) = c;
end

[a,lambda,exitflag] = linprog(f,Ac,bc,Aeq,beq,lb,[]);

A = helper_mat(a(2:end),N,M);

% Reconstruct 'source' signals S

switch( method )
    case 'inverse'
        S = pinv(A)*X;
    case 'mmse'
        noise_power = rel_noise_power * norm(X,'fro')^2;
        S = A'*inv(A*A' + noise_power*eye(N)) * X;
    case 'tls'
        [~,~,V] = svd([A X]);

        V_AX = V(1:M,M+1:end);
        V_XX = V(M+1:end,M+1:end);

        S = -V_AX / V_XX;
    otherwise
        error('source recovery method unsupported');
end
