function S = reconstruct_sources(A,X,method,tls_scaling)
% S = reconstruct_sources(A,X,method,tls_scaling)
%
% Estimates the sources from X and A(estimated).  

[N,M,K] = helper_get_network_size(X,A);

if ( nargin < 4 )
    tls_scaling = 1;
end

switch( method )
    case 'inverse'
        S = pinv(A)*X;
%     case 'mmse'
%         noise_power = rel_noise_power * norm(X,'fro')^2;
% %         S = A'*inv(A*A' + noise_power*eye(N)) * X; % don't compute
%         S = ( A' / (A*A' + noise_power*eye(N)) ) * X;
    case 'tls'
%         % naive approach
%         [~,~,V] = svd([A X]);
% 
%         V_AX = V(1:M,M+1:end);
%         V_XX = V(M+1:end,M+1:end);
% 
%         S = -V_AX / V_XX;
    
    % need to add scaling!
%     C = eye(M+K);
%     C = diag([10*ones(1,M) 0.01*ones(1,K)]);
    C = diag([tls_scaling*ones(1,M) 1*ones(1,K)]); % original: 10000 A, 1 X

    [Up,Ep,Vp] = svd(A*C(1:M,1:M));
    
    [U,E,V] = svd([A X]*C);
    
    U1 = U(:,1:M);
    U2 = U(:,M+1:end);
    
    E1 = E(1:M,1:M);
    E2 = E(M+1:end,M+1:end); % should be square? only need first entry.
    
    V11 = V(1:M,    1:M);
    V12 = V(1:M,    M+1:end);
    V21 = V(M+1:end,1:M);
    V22 = V(M+1:end,M+1:end);
    
    % check uniqueness criterion
    % want Ep(M,M) > E2(1,1)
    if ( Ep(M,M) <= E2(1,1) )
        warning('problem detected in TLS solution');
        display(sprintf('Ep(M,M) = %3.2e vs. E2(1,1) = %3.2e',Ep(M,M),E2(1,1)));
    end
    
%     S = - V12 * inv(V22); % compute implicitly
    S = tls_scaling * -V12 / V22;
    
%     AX = [A X]*C; % not needed
    
%     % debug
%     dAdX = U2 * E2 * [ V12' V22' ];
%     figure(10); clf;
%     surf(abs(dAdX));
%     ylabel(sprintf('N = %d',N));
%     xlabel(sprintf('M+K = %d+%d',M,K));
%     
%     avg_correction_A = norm(dAdX(:,1:M),'fro')^2 / (N*M)
%     avg_correction_X = norm(dAdX(:,(M+1):end),'fro')^2 / (N*K)
%     
%     1; % for breakpoint
%     % debug

    otherwise
        error('source recovery method unsupported');
end

