function [N,M,K] = helper_get_network_size(X,A_top)
% Does consistency checking and returns network size

[N,K] = size(X);

[N_A,M] = size(A_top);

if ( N_A ~= N )
    error('The sizes of A_top and X are inconsistent');
end
