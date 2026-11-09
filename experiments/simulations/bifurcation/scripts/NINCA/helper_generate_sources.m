function S = helper_generate_sources(M,K)
% S = helper_generate_sources(M,K) generates M source signals of length K
% default for K = 300

% default for K = 300
if ( nargin == 1 )
    K = 300;
end

S = zeros(M,K);

t = 1:K;

for m = 1:M
    pinch = 5 + 20*rand(1,1); % 5 + (0-20)
    offset = -K/4 + K/8*randn(1,1); % (-K/4) +/- (0-K/4)
    amplitude = 0.1 + 0.5*rand(1,1); % 0.1 + (0-0.5)
    
    S(m,:) = amplitude * exp( -((t+offset)/pinch).^2 );
end
