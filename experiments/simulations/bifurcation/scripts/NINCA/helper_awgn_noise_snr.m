function X_noisy = helper_awgn_noise_snr(X,snr_db)
% X_noisy = helper_awgn_noise_snr(X,snr_db)
% 
% Adds Additive White Guassian Noise (AWGN) to the data in X.
% This function assumes the X are noise-free, and adds noise to the
% observations.  The SNR is set by the parameter snr_db.

sig_power = norm(X,'fro')^2;

noise_power = sig_power / 10^(snr_db/10);

total_dim = size(X,1) * size(X,2);

N = randn(size(X)) ./ sqrt(total_dim) .* sqrt(noise_power);

% % debug
% actual_noise_power = norm(N,'fro')^2;
% keyboard;

X_noisy = X + N;
