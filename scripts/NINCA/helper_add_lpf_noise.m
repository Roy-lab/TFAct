function X_noisy = helper_add_lpf_noise(X,cutoff,rel_noise_power)
% X_noisy = helper_add_lpf_noise(X,cutoff,rel_noise_power)
% 
% Adds low pass filtered noise to the data in X.
% rel_noise_power determines the power of the noise in proportion to the
% power of the input.  E.g. rel_noise_power = 1.0 adds noise with the same
% total expected power of the input signal.  Specifically, X_noisy = X + N,
% where E{norm(N,'fro')^2} = rel_noise_power * norm(X,'fro')^2.
%
% cutoff in [0,1] denotes the digital cutoff frequency for the noise.

noise_power = rel_noise_power * norm(X,'fro')^2;

total_dim = size(X,1) * size(X,2);

N = randn(size(X)) ./ sqrt(total_dim) .* sqrt(noise_power); % white

[B,A] = ellip(6,0.5,20,cutoff); % ellip(order,Rp,Rs,cutoff)
% fvtool(B,A);

N = sqrt(1/cutoff) * filter(B,A,N,[],2); % filtered, should now be LP

% % debug
% actual_noise_power = norm(N,'fro')^2;
% keyboard;

X_noisy = X + N;
