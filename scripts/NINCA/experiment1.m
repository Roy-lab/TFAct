% Simulations for Oxford Bioinformatics Journal
%
% Experiment 1 data
%
% 'dataset2_bioinformatics1', 'dataset3_bioinformatics2',
% 'dataset4_bioinformatics3', 'dataset5_bioinformatics4'
%
% These are 4 different subnetworks of the larger networks/gene microarray
% data.  

%% Parameters

% choose algorithm
% algorithm = 'ncanec_quad';
% algorithm = 'ncanec_quad_nocvx';
% algorithm = 'ncanec_quad_nocvx_colwise'; % preferred
% algorithm = 'ncanec_lp_nocvx';
algorithm = 'ncanec_lp_nocvx_colwise'; % preferred
% algorithm = 'fastnca';
% algorithm = 'als_nocvx';
% algorithm = 'als';

% choose recovery method ( for necnca algorithms )
% method = 'inverse'; 
% method = 'mmse'; rel_noise_power = 0.05;
method = 'tls';
tls_scaling = 1e4; 

% truncate out last a-factor data
truncate = true; % use only first 3 time periods
% truncate = false; % use all 4

%% Run experiments

size_exp = zeros(4,3); % N, M, K for each of 4 experiements

M = 40;

if ( truncate ) 
    K = 56;
else
    K = 69;
end

size_exp(1,:) = [1247; M; K];
size_exp(2,:) = [943; M; K];
size_exp(3,:) = [1144; M; K];
size_exp(4,:) = [921; M; K];

A = zeros(4,max(size_exp(:,1)),M); % A is 4 x N_max x M
S = zeros(4,M,K); 

for exp = 1:4
    
    exp
    
    % load data, provides X, A_topology
    switch (exp)
        case 1
            load('dataset2_bioinformatics1');
        case 2
            load('dataset3_bioinformatics2');
        case 3
            load('dataset4_bioinformatics3');
        case 4
            load('dataset5_bioinformatics4');
    end
    
    % truncate
    X = X(:,1:K);
    
    % run algorithm to recover A,S
    if ( ~exist('rel_noise_power','var') )
        rel_noise_power = 0;
    end
    
    tic;
    
    switch ( algorithm )
        case 'als'
            [Ai,Si] = als_sub( X, A_topology );
        case 'als_nocvx'
            [Ai,Si] = als_sub_nocvx( X, A_topology );
        case 'fastnca'
            [Ai,Si] = FastNCA( X, A_topology );
%     case 'ncanec_quad'
%         [Ai,Si] = ncanec_quad_sub( X, A_topology, rel_noise_power, method );
%     case 'ncanec_quad_nocvx'
%         [Ai,Si] = ncanec_quad_sub_nocvx( X, A_topology, rel_noise_power, method );
    case 'ncanec_quad_nocvx_colwise'
        [Ai,Si] = ncanec_quad_sub_nocvx_colwise( X, A_topology, method, tls_scaling );
%     case 'ncanec_lp_nocvx'
%         [Ai,Si] = ncanec_lp_sub_nocvx( X, A_topology, rel_noise_power, method );
    case 'ncanec_lp_nocvx_colwise'
        [Ai,Si] = ncanec_lp_sub_nocvx_colwise( X, A_topology, method, tls_scaling );
        otherwise
            error('Unrecognized algorithm requested');
    end
    
    toc
    display(sprintf('runtime on exp %d using algo %s',exp,algorithm));
    
    % save into struct
    A(exp,1:size_exp(exp,1),1:M) = Ai;
    S(exp,:,:) = Si;
    
end


%% Plot figure

% TFAs of interest ( 0 => not present )
tfa_index = [  2  9 10 15 16 22 31  0 34 35 36; ... % exp 1
               2  0  8 18 19 22 30  0 33 34 35; ... % exp 2
               2  0  7 16 17 20 30  0 32 33 34; ... % exp 3
               1  0  6 18 19 22 30 32 35 36 37 ];   % exp 4
           

% draw plot
figure; clf;

spi = 0;
for tfa = 1:11    
    % read tfa into buffer ( convert 1x1xK -> 1xK )
    S_buffer = zeros(4,K); % exp x K
    for exp = 1:4
        for k = 1:K
            if ( tfa_index(exp,tfa) ~= 0 )
                S_buffer(exp,k) = S(exp,tfa_index(exp,tfa),k);
            end
        end
    end
    
    if ( truncate )
        time_periods = 3;
    else
        time_periods = 4;
    end
    
    for tp = 1:time_periods 
        spi = spi + 1;
        
        subplot(11,time_periods,spi);
        hold on;
        
        switch (tp)
            case 1 % elutriation
                t = 0:30:390;
                if ( tfa_index(1,tfa) ~= 0 )
                    plot( t, S_buffer(1,1:14), 'b-');
                end
                if ( tfa_index(2,tfa) ~= 0 )
                    plot( t, S_buffer(2,1:14), 'b:');
                end
                if ( tfa_index(3,tfa) ~= 0 )
                    plot( t, S_buffer(3,1:14), 'b-.');
                end
                if ( tfa_index(4,tfa) ~= 0 )
                    plot( t, S_buffer(4,1:14), 'b--');
                end
                
                % set axis nicely
                v = axis;
                axis([0 400 v(3) v(4)]);
                
                if ( tfa == 1 )
                    title('elutriation');
                end
            case 2 % \alpha-factor
                t = 0:7:119;
                if ( tfa_index(1,tfa) ~= 0 )
                    plot( t, S_buffer(1,15:32), 'b-');
                end
                if ( tfa_index(2,tfa) ~= 0 )
                    plot( t, S_buffer(2,15:32), 'b:');
                end
                if ( tfa_index(3,tfa) ~= 0 )
                    plot( t, S_buffer(3,15:32), 'b-.');
                end
                if ( tfa_index(4,tfa) ~= 0 )
                    plot( t, S_buffer(4,15:32), 'b--');
                end
                
                % set axis nicely
                v = axis;
                axis([0 120 v(3) v(4)]);
                
                if ( tfa == 1 )
                    title('\alpha-factor');
                end
            case 3 % cdc-15
                t = [10 30 50 70:10:250 270 290];
                if ( tfa_index(1,tfa) ~= 0 )
                    plot( t, S_buffer(1,33:56), 'b-');
                end
                if ( tfa_index(2,tfa) ~= 0 )
                    plot( t, S_buffer(2,33:56), 'b:');
                end
                if ( tfa_index(3,tfa) ~= 0 )
                    plot( t, S_buffer(3,33:56), 'b-.');
                end
                if ( tfa_index(4,tfa) ~= 0 )
                    plot( t, S_buffer(4,33:56), 'b--');
                end
                
                % set axis nicely
                v = axis;
                axis([0 300 v(3) v(4)]);
                
                if ( tfa == 1 )
                    title('cdc-15');
                end
            case 4 % \alpha-factor (mutant)
                t = 0:15:180;
                if ( tfa_index(1,tfa) ~= 0 )
                    plot( t, S_buffer(1,57:69), 'b-');
                end
                if ( tfa_index(2,tfa) ~= 0 )
                    plot( t, S_buffer(2,57:69), 'b-');
                end
                if ( tfa_index(3,tfa) ~= 0 )
                    plot( t, S_buffer(3,57:69), 'b-');
                end
                if ( tfa_index(4,tfa) ~= 0 )
                    plot( t, S_buffer(4,57:69), 'b-');
                end
                
                % set axis nicely
                v = axis;
                axis([0 180 v(3) v(4)]);
                
                if ( tfa == 1 )
                    title('\alpha-factor (mutant)');
                end
        end
        
        
    end
end

% viewer for individual sources
exp1data_viewer('K',K,'S',S,'tfa_index',tfa_index)


            