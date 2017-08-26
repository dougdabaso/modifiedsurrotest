% dist_locvsglob.m
%
% Jun Xiao, Pierre Borgnat & Patrick Flandrin 
% 09/2007
%
% computes the distances between local TFR and frequency marginal mean of that TFR
%
% - input:   - S : TFR (positive frequency only)
%            - t : time vector  
%            - [opt.] opt2 : option of distance (1<=opt2<=8)
%               1-Lq, 2-Kolmogorov, 3-Kullback-Leibler, 4-Jensen Divergence
%               5-log-spectral deviation, 6-Itakura-Saito, 
%               7-maison1: Ku*(1+abs(log(sum(S1)/sum(S2)))), 
%               8-maison2: computation of distance maison Ku*(1+LSD) 
%               (default 8)
%            - [opt.] a & b : min frequency & max frequency (0=<a<b=<0.5) 
%               (default: a=0, b=0.5) 
%            - [opt.] lambda : parameter of ad hoc distance 
%               (default: lambda=1)
%               9 - diffusion distance
%              10 - Generalized Matusita distance (for r = 2)
%              11 - Symmetrized itakuro
%
% - output:  - dS : the distances between marginals and multitapers of spectrogram
%
% Usage: dS = dist_locvsglob(S,t,opt2,a,b) ;


% [was dist_MTss2.m]

function dS = dist_locvsglob(S,t,opt2,a,b,lambda) ;

%a = 0;
%b = 0.5;
lambda = 1;

ymargS=[]; ytapS=[]; dS=[]; 
q = 2;  

[Nfftover2,Nx]=size(S);
B = max(round(a*Nfftover2*2),1):round(b*Nfftover2*2);  % frequency band 


ymargS = mean(S,2) ; % mean of marginals

if opt2 == 1
    
    ymargSN = ymargS/sum(ymargS) ;  % normalization 1
    
else
    
    ymargSN = abs(ymargS)/sum(abs(ymargS)) ; % normalization 2
    
end

for n = 1:length(t)

    ytapS = S(:,n) ; % mean of tapers 
    
    if opt2 == 1 
        
        ytapSN = ytapS/sum(ytapS) ; %  normalization 1  
        dS(n) = (sum(abs(ymargSN(B,:)-ytapSN(B,:)).^q)).^(1/q); % computation of distance Lq
        
    elseif opt2 == 2
        
        ytapSN = abs(ytapS)/sum(abs(ytapS)) ; % normalization 2
        dS(n) = sum(abs(ymargSN(B,:)-ytapSN(B,:)));  % computation of distance Ko 
        
        
    elseif opt2 == 3
        
        ytapSN = abs(ytapS)/sum(abs(ytapS)) ; % normalization 2       
        %dS(n) = sum((ytapSN(B,:)-ymargSN(B,:)).*log(ytapSN(B,:)./ymargSN(B,:))) ; % computation of distance Ku
        dS(n) = sum((ytapSN(B,:)).*log(ytapSN(B,:)./ymargSN(B,:))) ; % computation of distance Ku
        
        
    elseif opt2 == 4
        
        ytapSN = abs(ytapS)/sum(abs(ytapS)) ; % normalization 2
        J1 = (ytapSN(B,:)+ymargSN(B,:))/2 ; % computation of distance  Je
        J2 = ytapSN(B,:) ;
        dS(n) = Rp(sqrt(J1.*J2))-(Rp(J1)+Rp(J2))/2 ;   
        
    elseif opt2 == 5
        
        dS(n) = sum(abs(log(ytapS(B,:)./ymargS(B,:)))) ; % computation of distance q (Log spectral deviation) 
        
    elseif opt2 == 6
        
        dS(n) = sum(abs(ytapS(B,:)./ymargS(B,:) - log(ytapS(B,:)./ymargS(B,:)) - 1)) ; % computation of distance Is
        


    elseif opt2 == 7
        
        ytapSN = abs(ytapS)/sum(abs(ytapS)) ; % normalization 2       
        dS(n) = sum((ytapSN(B,:)-ymargSN(B,:)).*log(ytapSN(B,:)./ymargSN(B,:))).*(1 + lambda*abs(log(sum(ytapS(B,:))/sum(ymargS(B,:))))) ; % computation of distance maison Ku*(1+abs(log(sum(S1)/sum(S2))))
        
elseif opt2 == 8
        
        ytapSN = abs(ytapS)/sum(abs(ytapS)) ; % normalization 2       
        dS(n) = sum((ytapSN(B,:)-ymargSN(B,:)).*log(ytapSN(B,:)./ymargSN(B,:))).*(1 + lambda*sum(abs(log(ytapS(B,:)./ymargS(B,:)))) ) ; % computation of distance maison Ku*(1+Lq)     

 elseif opt2 == 9
     
 dS(n) = diffusion_distance(ytapS(B,:),ymargS(B,:));
  
  elseif opt2 == 10
       

         ytapSN = abs(ytapS)/sum(abs(ytapS)) ; % normalization 2
        dS(n) =  sum((sqrt(ymargSN(B,:))-sqrt(ytapSN(B,:))).^2);  % computation of distance Ko 
 
  elseif opt2 == 11
     
 
%symmetrized itakuro
        %dS(n) = (sum(abs(ytapS(B,:)./ymargS(B,:) - log(ytapS(B,:)./ymargS(B,:)) - 1)) + sum(abs(ymargS(B,:)./ytapS(B,:) - log(ymargS(B,:)./ytapS(B,:)) - 1)))/2;
        dS(n) = sum(abs(ytapS(B,:)./ymargS(B,:) + ymargS(B,:)./ytapS(B,:) - 2)) ; % computation of distance Is
%         size(ytapS)    
        
    end    

end

% plot(dS)
% hold on



return;

