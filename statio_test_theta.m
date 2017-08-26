% statio_test_theta.m
%
% Jun Xiao, Pierre Borgnat & Patrick Flandrin 
% 09/2007
%
% Computation of the Statistics Theta for Stationarity Test
%
% Input:  - tfr: matrix of TFRepresentation, positive frequency only
%         - tt2: time vector used
%         - [opt] opt_dist: distance used in dist_locvsglob
%         - [opt] a,b: frequency parameters of dist_locvsglob
%       
% Output:  - theta: statistics used for Stationarity Test
%
% Need: dist_locvsglob
%
% Usage: [theta,Cn_dist,Cn_mean]=statio_test_theta(tfr,tt2,opt_dist,a,b);


function [theta,Cn_dist,Cn_mean]=statio_test_theta(tfr,tt2,opt_dist,a,b);
%  length(tt2)

[l,c]=size(tfr);

tt2 = 1:c;
%a = 0;
%b = 0.5;
lambda=1;

Cn_dist = dist_locvsglob(tfr,tt2,opt_dist,a,b,lambda) ; 
Cn_mean = mean(Cn_dist); 

theta = sum((Cn_dist - Cn_mean).^2)/(length(Cn_dist)-1); 


return;
