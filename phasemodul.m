% phasemodul.m
%
% Jun Xiao, Pierre Borgnat & Patrick Flandrin 
% 09/2007
%
% Replaces the phase of determined signal by a random one (uniform in [0 2*pi]).
%
% - input  : - x : determined signal
%            - [opt] Nx : length of the FFT transform 
%              (default:  Nx=length(x))
%
% - output : - z : signal with the same amplitude of signal x but an random phase
%
% Usage: z = phasemodul(x,Nx);

function z = phasemodul(x,Nx);
 
% if nargin==1,
%     Nx=length(x);
% end;

Nfft = Nx ;
%y = [] ; 
A = zeros(Nfft,1) ; 
%phase0 = []; z=[];  

y = fft(x,Nfft);
A(1:ceil(Nfft/2)) = abs(y(1:ceil(Nfft/2)));

phase0 = rand(Nfft,1).*2.*pi;
z = real(ifft(2*A.*exp(i.*phase0)));

if size(z)~=size(x),
    z=z.';
end;

return;
    