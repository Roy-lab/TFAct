function X_noisy = helper_add_space_corr_noise(X,corr_fac,rel_noise_power)
% X_noisy = helper_add_lpf_noise(X,corr_fac,rel_noise_power)
% 
% Adds spatially correlated noise to X.  Correlation controlled by corr_fac
% in [0,1].
% 
% rel_noise_power determines the power of the noise in proportion to the
% power of the input.  E.g. rel_noise_power = 1.0 adds noise with the same
% total expected power of the input signal.  Specifically, X_noisy = X + N,
% where E{norm(N,'fro')^2} = rel_noise_power * norm(X,'fro')^2.


noise_power = rel_noise_power * norm(X,'fro')^2;

total_dim = size(X,1) * size(X,2);

N = randn(size(X)) ./ sqrt(total_dim) .* sqrt(noise_power); % white

corr_vec = corr_fac.^(0:size(X,1)-1);

corr_mat = toeplitz(corr_vec);

[V,D] = eig(corr_mat);

% form correlating filter
C = V*sqrt(D);
        
% color noise
N = C * N;

% % debug
% actual_noise_power = norm(N,'fro')^2;

X_noisy = X + N;
