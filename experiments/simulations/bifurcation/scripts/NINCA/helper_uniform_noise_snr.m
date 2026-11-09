function X_noisy = helper_uniform_noise_snr(X,snr_db,positive_only)
% X_noisy = helper_uniform_noise_snr(X,snr_db)
% 
% Adds uniformly distributed noise to the data in X.  
% This function assumes the X are noise-free, and adds noise to the
% observations.  The SNR is set by the parameter snr_db.  The parameter
% positive_only controls whether the noise is stricly positive or not. 
% Default is false (i.e. noise is positive and negative).

if ( nargin < 3 )
    positive_only = false;
end

sig_power = norm(X,'fro')^2;

noise_power = sig_power / 10^(snr_db/10);

total_dim = size(X,1) * size(X,2);

if ( positive_only )
    N = rand(size(X))       ./ sqrt(total_dim) .* sqrt(3*noise_power);
else
    N = (rand(size(X))-0.5) ./ sqrt(total_dim) .* sqrt(12*noise_power);
end

% % debug
% noise_power
% actual_noise_power = norm(N,'fro')^2
% keyboard;

X_noisy = X + N;
