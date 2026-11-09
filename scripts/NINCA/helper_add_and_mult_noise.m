function X_noisy = helper_add_and_mult_noise(X_true,noise_scale)
% X_noisy = helper_add_and_mult_noise(X_true,noise_scale)
% 
% Implements noise model from Kossendov and Ochs, "Matrix Factorization
% Methods Applied in Microarray Data Analysis".

if ( nargin < 2 )
    noise_scale = 1; 
end

[N,K] = size(X_true);

sig_power = norm(X_true,'fro')^2/numel(X_true);

% Parameters from Kossenkov paper
N_trials = 4;
sig_a = noise_scale*0.16*sqrt(max(max(abs(X_true)))); % additive noise std. dev.
sig_b = noise_scale*0.29*sqrt(sig_power); % multiplicative noise std. dev.


D = zeros(size(X_true));

for trial = 1:N_trials
    for n = 1:N
        for k = 1:K
            D(n,k) = D(n,k) + sig_a*randn(1) + X_true(n,k)*exp(sig_b*randn(1));
        end
    end
end

X_noisy = D / N_trials;
