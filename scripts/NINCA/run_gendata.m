% Run a dataset using a particular algorithm and noise setting.

clear;

disp(['Running run_gendata.m @ ' datestr(now) '..']);

%  PARAMETERS ----------------------------------------------------------

% Choose algorithm
% algorithm = 'als'
% algorithm = 'als_nocvx'
% algorithm = 'fastnca'
% algorithm = 'ncanec_quad'
% algorithm = 'ncanec_quad_nocvx'
% algorithm = 'ncanec_quad_nocvx_colwise'
% algorithm = 'ncanec_lp_nocvx'
algorithm = 'ncanec_lp_nocvx_colwise'

% Choose noise type
noise_profile = 'awgn'; snr_db = 25;
% noise_profile = 'quantization'; 


% RUN ------------------------------------------------------------------

% Load data
load( 'A_true' ); % A_true

A_connectivity = double(A_true > 0);

N = 7; M = 3;

% load( 'dataset_S_highcorr1' ); % S
load( 'dataset_S_highcorr2' ); % S
% load( 'dataset_S_lowcorr1' ); % S
% load( 'dataset_S_mixedcorr1' ); % S
S_true = S;

X = A_true*S_true;

% Add noise
switch( noise_profile )
    case 'awgn'
        X = helper_awgn_noise_snr(X,snr_db);
    case 'quantization'
        error('The quantization noise option is not implemented yet');
    otherwise
        error('Unrecognized noise profile');
end

% Run Algorithm
switch ( algorithm )
    case 'als'
        [A,S] = als_sub( X, A_connectivity );
    case 'als_nocvx'
        [A,S] = als_sub_nocvx( X, A_connectivity );
    case 'fastnca'
        [A,S] = FastNCA( X, A_connectivity );
    case 'ncanec_quad'
        [A,S] = ncanec_quad_sub( X, A_connectivity );
    case 'ncanec_quad_nocvx'
        [A,S] = ncanec_quad_sub_nocvx( X, A_connectivity );
    case 'ncanec_quad_nocvx_colwise'
        [A,S] = ncanec_quad_sub_nocvx_colwise( X, A_connectivity );
    case 'ncanec_lp_nocvx'
        [A,S] = ncanec_lp_sub_nocvx( X, A_connectivity );
    case 'ncanec_lp_nocvx_colwise'
        [A,S] = ncanec_lp_sub_nocvx_colwise( X, A_connectivity );
    otherwise
        error('Unrecognized algorithm requested');
end

% Overall Reconstruction Error
reconstruction_error = norm(X-A*S,'fro')

% Report MSE of A recovery, if we have that information
if ( exist('A_true','var') )
    
    % Scale and report connectivity matrix A
    A_true_scaled = zeros(size(A_true));
    A_scaled = zeros(size(A));

    for c = 1:M
        A_true_scaled(:,c) = A_true(:,c) ./ sum(A_true(:,c));
        A_scaled(:,c) = A(:,c) ./ sum(A(:,c));
    end

    MSE_A = norm( A_true_scaled - A_scaled, 2 )^2 / (N*M)

end

% Plot source recovery
figure(5); clf;
set(gcf,'Name',sprintf('source recovery for %s',algorithm));

good_axis = [0 300 -0.25 0.75];

% original signals
subplot(3,2,1); plot(S_true(1,:)); axis(good_axis);
subplot(3,2,3); plot(S_true(2,:)); axis(good_axis);
subplot(3,2,5); plot(S_true(3,:)); axis(good_axis);

% recovered signals

subplot(3,2,2); plot(S(1,:)); axis(good_axis);
subplot(3,2,4); plot(S(2,:)); axis(good_axis);
subplot(3,2,6); plot(S(3,:)); axis(good_axis);

