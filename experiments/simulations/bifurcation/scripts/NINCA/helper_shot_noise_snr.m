function X_noisy = helper_shot_noise_snr(X,prob_of_occurence,snr_db)
% X_noisy = helper_shot_noise_snr(X,prob_of_occurence,snr_db)
% 
% Adds Shot noise as described in helper_add_shot_noise.m .
% This function assumes the X are noise-free, and adds noise to the
% observations.  The SNR is set by the parameter snr_db.

sig_power = norm(X,'fro')^2;

noise_power = sig_power / 10^(snr_db/10);

n = size(X,1);
k = size(X,2);

total_dim = n*k;

one_total_dim = 1/total_dim;
one_prob_of_occurence = 1/prob_of_occurence;

N = zeros(size(X));
for i = 1:n
    for j = 1:k
        if ( rand(1) < prob_of_occurence ) % insert noise on this sample
            N(i,j) = sqrt( noise_power * one_total_dim * one_prob_of_occurence ) * randn(1);
        end
    end
end

% % debug
% noise_power
% actual_noise_power = norm(N,'fro')^2
% keyboard;

X_noisy = X + N;
