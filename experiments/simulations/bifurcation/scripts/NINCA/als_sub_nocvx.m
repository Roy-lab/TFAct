function [A,S] = als_sub_nocvx( X, A_top )
% [A,S] = als_sub( X, A_top ) computes the recovered connectivity matrix A
% and the recovered source signals S given observations X.  N and M are
% the number of genes and TFs, respecitively, with the convention that A is
% an N x M matrix.  

% parameters 
MAX_ITERATIONS = 100; % if convergence fails, don't run forever

% find N, M, K; check consistency
[N,M,K] = helper_get_network_size(X,A_top);


% find zero locations in A from A_top
c = 0.001; % "zero" threshold
I = find( A_top  < c );
J = find( A_top >= c );


% run ALS algorithm

% ALS
% initial A: just use normalized A_top as a starting point
% start ALS by using this A and finding an S

A = A_top / max(max(A_top)); % make largest entry == 1
S = zeros(M,K);

error = norm( X - A*S, 'fro' );
new_error = 100*error; % initialize to something much larger than error

it = 0;
while ( (abs(error-new_error) > 0.001*new_error) && it < MAX_ITERATIONS )
    
    it = it + 1;
    
    error = new_error;
    
    % fix A, update S
%     cvx_begin quiet
%         variable S(M,K);
% 
%         minimize( norm( X - A*S, 'fro' ) );
%     cvx_end
    H = kron(sparse(eye(K)),A'*A);
    f = -vec(A'*X);
    
%     H_S_percent0 = sum(sum(H==0))/prod(size(H)) % debug
%     H = sparse(H);
    
    s = quadprog(H,f,[],[]);
    
    S = helper_mat(s,M,K);
    
    % fix S, update A
%     cvx_begin quiet
%         variable A(N,M);
% 
%         minimize( norm( X - A*S, 'fro' ) );
% 
%         subject to
%             % zero locations
%             for i = I
%                 A(I) == 0;
%             end
%     cvx_end
    H = kron(S*S',sparse(eye(N)));
    f = -vec(X*S');
    
    % A(I) == 0
    Aeq = zeros(length(I),N*M);
    beq = zeros(length(I),1);
    A_0 = zeros(size(A));
    for i_index = 1:length(I)
        A_0(I(i_index)) = 1;
        Aeq(i_index,:) = vec(A_0)';
        
        A_0(I(i_index)) = 0; % reset entry to zero
    end
    
%     H_A_percent0 = sum(sum(H==0))/prod(size(H)) % debug
%     H = sparse(H);
%     Aeq = sparse(Aeq);
    
    a = quadprog(H,f,[],[],Aeq,beq);
    
    A = helper_mat(a,N,M); % missing before??!? (3/8/11)
    
    new_error = norm(X-A*S,'fro');
end

% check if we hit max iterations
if ( it == MAX_ITERATIONS )
    display('WARNING: ALS did not converge');
else
    display(sprintf('ALS algorithm completed in %d iterations',it));
end


