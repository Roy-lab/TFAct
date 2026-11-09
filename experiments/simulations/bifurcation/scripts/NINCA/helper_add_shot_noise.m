function X_noisy = helper_add_shot_noise(X,prob_of_occurence,rel_noise_power)
% X_noisy = helper_add_shot_noise(X,rel_noise_power)
% 
% Adds shot noise to the data in X.
% rel_noise_power determines the power of the noise in proportion to the
% power of the input.  E.g. rel_noise_power = 1.0 adds noise with the same
% total expected power of the input signal.  Specifically, X_noisy = X + N,
% where E{norm(N,'fro')^2} = rel_noise_power * norm(X,'fro')^2.
%
% For the purposes of this function, we'll add shot noise like this:
% noise = { sqrt( (1/prob_of_occurence)*noise_power )*randn(1) w/ prob = prob_of_occurence;
%         { 0                                                  w/ prob = 1 - prob_of_occurence

noise_power = rel_noise_power * norm(X,'fro')^2;

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
