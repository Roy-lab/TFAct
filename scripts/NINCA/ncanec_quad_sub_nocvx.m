function [A,S] = ncanec_quad_sub_nocvx( X, A_top, rel_noise_power, method )
% [A,S] = ncanec_quad_sub( X, A_top ) computes the recovered connectivity matrix A
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

% cvx_begin
%     variable A(N,M);
%     
%     minimize( norm(U_0'*A,'fro') );
%     
%     subject to
%         % zeros
%         for i = I
%             A(i) == 0;
%         end
%         
%         % non-negativity
%         for j = J
%             A(j) >= c;
%         end
%         
%         % avoid nontrivial solution
%         for m = 1:M
%             sum(A(:,m)) == 1;
%         end
%         
% cvx_end
% 
% A = full(A); % if A is sparse, convert it to a full matrix for pinv

H = sparse( kron(sparse(eye(M)),U_0*U_0') ); 

% A(I) == 0 
Aeq0 = zeros(length(I),N*M);
beq0 = zeros(length(I),1);
A_0 = zeros(size(A_top));
for i_index = 1:length(I)
    A_0(I(i_index)) = 1;
    Aeq0(i_index,:) = vec(A_0)';

    A_0(I(i_index)) = 0; % reset entry to zero
end

% sum(cols(A)) == 1
Aeq1 = zeros(M,N*M);
beq1 = ones(M,1);

for m = 1:M
    A_0 = zeros(size(A_top));
    A_0(:,m) = 1;
    Aeq1(m,:) = vec(A_0)';
end

% build Aeq, beq
Aeq = [ Aeq0; Aeq1 ]; 
beq = [ beq0; beq1 ];

Aeq = sparse(Aeq);

% A(J) >= c
lb = zeros(N*M,1);
[J_r,J_c] = ind2sub(size(A_top),J);
for j_index = 1:length(J)
    lb( J_r(j_index) + N*(J_c(j_index)-1) ) = c;
end

[a,fval,exitflag] = quadprog(H,[],[],[],Aeq,beq,lb,[],vec(A_top));

A = helper_mat(a,N,M);

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


