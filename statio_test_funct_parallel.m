%%%%% Modified Stationarity Test Function, version 1.12 %%%%%%%%%%%%%%%%
%                                                                      %
% This algorithm tests stationarity using TF frequency representation, %
% stationarity surrogates, distance-based and weighted techniques and  %
%  using the time-marginal-weighting procedure described in:           %
%                                                                      %
% D. B. de Souza, J. Chanussot, A. C. Favre and P. Borgnat,"A modified %
% time-frequency method for testing wide-sense stationarity," 2012     %
% IEEE International Conference on Acoustics, Speech and Signal        %
% Processing (ICASSP), Kyoto, 2012, pp. 3409-3412.                     %
%                                                                      %
% The original algorithm belongs to P. Borgnat, J. Xiao and P.Flandrin %
%                                                                      %
% References: see perso.ens-lyon.fr/pierre.borgnat                     %
%                                                                      %
% 1) "Testing Stationarity with Time-Frequency Surrogates",            %
% J. Xiao,P. Borgnat,P. Flandrin, 15th European Signal Processing      %
% Conference EUSIPCO-2007,  Poznan (Pologne), 3-7 septembre 2007       %
% 2)  "Sur un test temps-fr?quence de stationnarit?",J.Xiao,P.Borgnat, %
% P. Flandrin, 21e Colloque sur le Traitement du Signal et des Images. %
% GRETSI-2007, Troyes (France), 11-14 septembre 2007.                  %
%                                                                      %
%                                                                      %
% Douglas David Baptista de Souza - 2014                               %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Input: x_array:    signal or array of signals to test
%        Nh0:        freq-time tradeoff resolution parameter (default Nh0 = 0.1)
%        fa:         parameter for computing theta dist. (default fa = 0)
%        fb:         parameter for computing theta dist. (default fa = 0.5)
%        opt_dist:   distance used by the method (default opt_dist = 6)
%        JJ:         Number of surrogates used to build null distribution (default JJ = 75)


%
% Output (inside the struct test_results_struct):

%           theta1            Test Statistics for original signal
%           theta0            Test Statistics for each surrogate (vector of size opttest.JJ)
%           gam_hat           Estimation of [alpha,beta] of the Gamma distribution of the {Theta_0}
%           gam_thresh        Threshold of the test
%           res_test          Result of the test: 0=stationary, 1=non-stationary
%           INS               Index of Non-stationarity
%           theta0_weigthed   Vector theta0 after applying the weighting method
%           gam_hat2          Estimation of [alpha,beta] of the Gamma distribution of the theta0_weigthed
%           gam_thresh2       Threshold of the test after weighting
%           res_test2         esult of the test after weighting
%           INS2              Index of Non-stationarity after weighting
%
%
%
%
% Need: phasemodul dist_locvsglob statio_test_theta tfrsp_h tfrsp_hm
% mean_hmt5 obtaining_gamthresh


%% Distances

%% Frequency-based
% 6 - Itakura-Saito
% 5 - log-spectral deviation
% 9 - Diffusion

%% Probability-based
% 3 - Kullback-Leibler
% 2 - Kolmogorov
% 10 - Generalized Matusita distance (for r = 2)

%% Mixed domains
% 8 - Combined
% 11 - Symmetrized itakuro



function test_results_struct = statio_test_funct_parallel(x_array,Nh0,fa,fb,opt_dist,JJ);


% Which average will be taken when computing the multitaper spectrogram
% default (mean) opt1=1
opt1=1;

% Setting the false alarm rate
fa_rate=0.05;

Nx = size(x_array,2);
number_real = size(x_array,1);

% For all signals in the array, set mean to zero and std to one
%parfor i = 1:size(x_array,1)
for i = 1:size(x_array,1)
    x_array(i,:) =  x_array(i,:) - mean(x_array(i,:));
    x_array(i,:) = x_array(i,:)/std(x_array(i,:));
end

% Setting parameters for the time-frequency representation
Nh=round(Nx*Nh0/2)*2-1;
Nfft=2^ceil(log2(Nh));
dt=fix((Nh+1)/8);
sides=(Nh+1)/2;
tt=sides:dt:Nx-sides;
ttred=tt;
Mh=5;
tm=5;

% Computation of multitaper spectrogram

%parfor i = 1:number_real
for i = 1:number_real
    MSp = tfrsp_hm(hilbert(x_array(i,:)).',tt,Nfft,Nh,Mh,tm);
    tfrx = mean_hmt5(MSp,opt1);
    tfr=tfrx(1:Nfft/2,:);
    time_marginal = sum(tfr,1)./size(tfr,1);
    norm_time_marginal(i,:) = time_marginal/max(time_marginal);
    [theta1(i),Cn_dist,Cn_mean]=statio_test_theta(tfr,ttred,opt_dist,fa,fb);
    Cn_dist_weigthed = Cn_dist.*norm_time_marginal(i,:);
    new_theta1(i) = sum((Cn_dist_weigthed - mean(Cn_dist_weigthed)).^2)/(size(Cn_dist_weigthed,2)-1);
end


% obtaining the test statistic for the surrogate set

%parfor jj=1:JJ*number_real
i_real = 1;
for jj=1:JJ*number_real
    % Replacing the phase of the signal by a random one
    % z = phasemodul(x_array(floor(JJ/jj)+1,:),Nx);
    % array_z(jj,:)=z;
    
    z = phasemodul(x_array(i_real,:),Nx);
    MSp = tfrsp_hm(hilbert(z),tt,Nfft,Nh,Mh,tm);
    tfrz = mean(MSp,3);
    tfr=tfrz(1:Nfft/2,:);
    [theta0(jj),Cn_distz,Cn_meanz]=statio_test_theta(tfr,ttred,opt_dist,fa,fb);
    Cn_dist_surrogates_weigthed = Cn_distz.*norm_time_marginal(i_real,:);
    theta0_weigthed(jj) = sum((Cn_dist_surrogates_weigthed - mean(Cn_dist_surrogates_weigthed)).^2)/(size(Cn_dist_surrogates_weigthed,2)-1);

    if jj == JJ
    i_real = i_real + 1;    
    end
    
end
clear norm_time_marginal


for k = 1:number_real
    up_bound = k*JJ;
    lo_bound = (k-1)*(JJ)+1;
    
    gam_hat = gamfit(theta0(lo_bound:up_bound),'mle');
    gam_hat2 = gamfit(theta0_weigthed(lo_bound:up_bound),'mle');
    gam_thresh=gaminv(1-fa_rate,gam_hat(1),gam_hat(2));
    gam_thresh2=gaminv(1-fa_rate,gam_hat2(1),gam_hat2(2));
    % Checking the test outcome, before weighting
    res_test = (theta1(k)>gam_thresh);
    % Checking the test outcome, after weighting
    res_test2 = (new_theta1(k)>gam_thresh2);
    % Gathering the INS value, before weighting
    INS = sqrt(theta1(k)/mean(theta0(lo_bound:up_bound)));
    % Gathering the INS value, after weighting
    INS2 = sqrt(new_theta1(k)/mean(theta0_weigthed(lo_bound:up_bound)));
    
    
    
    test_results_struct.(sprintf('real_%d',k)) = struct('theta0',theta0(lo_bound:up_bound),...
        'theta0_weigthed',theta0_weigthed(lo_bound:up_bound),...
        'theta1',theta1(k),...
        'new_theta1',new_theta1(k),...
        'gam_hat',gam_hat,...
        'gam_hat2',gam_hat2,...
        'gam_thresh',gam_thresh,...
        'gam_thresh2',gam_thresh2,...
        'INS',INS,...
        'INS2',INS2,...
        'res_test',res_test,...
        'res_test2',res_test2)
    
end

clear theta0
clear theta0_weigthed
clear theta1
clear new_theta1
clear gam_hat
clear gam_hat2
clear gam_thresh
clear gam_thresh2
clear INS
clear INS2
clear res_test
clear res_test2





