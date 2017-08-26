% tfrsp_hm.m
%
% P. Flandrin & J. Xiao, 2005
%
% computes Hermite multitaper (reassigned) spectrograms
%
% input  - x : signal
%        - t : time instants of analysis
%        - Nfft : number of frequency bins for FFT
%        - Nh : number of points for Hermite tapers (must be odd) 
%        - M : maximum order 
%        - tm : half time support (>= 6 recommended)
%
% output - S : spectrograms (Nfft x length(t) x M)
%        - RS : reassigned spectrograms (Nfft x length(t) x M)
%        - hat : reassigned vector fields (Nfft x length(t) x M)
%        - tt : Hermite support (1 x Nh)
%
% calls  - tfrrsp_h.m
%        - hermf.m

function S = tfrsp_hm(x,t,Nfft,Nh,M,tm) ;

if size(x,1) == 1
   x = x.';
end

[h,Dh,tt] = hermf(Nh,M,tm) ;

S = zeros(Nfft,length(t),M) ; 

for k = 1:M    
    spt = tfrsp_h(x,t,Nfft,h(k,:)',Dh(k,:)') ;
    
    S(:,:,k) = spt ;
    
end