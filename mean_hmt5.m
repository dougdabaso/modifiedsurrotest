% mean_hmt5.m
%
% P. Flandrin, Nov. 18, 2005
%
% computes various means of multitaper spectrograms
%
% input  - S : spectrograms (Nfft x Nx x M)
%        - opt : option for type of means, 1 = arithmetic (default),
%          2 = geometric, 3 = min, 4 = median
%
% output - Sm mean spectrogram (Nfft x Nx)

function Sm = mean_hmt5(S,opt) ;

if opt == 1
    
    Sm = mean(S,3) ;
    
elseif opt == 2
    
    prec = 1e-14 ;
    Sm = exp(mean(log(max(S,prec)),3)) ;
    
elseif opt == 3
    
    Sm = min(S,[ ],3) ;

elseif opt == 4
    
    Sm = median(S,3) ;
    
end


