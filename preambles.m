
%%
CS = 20e6; % channel spacing
NSD = 48; %number of data subcarriers
NSP = 4; % number of pilot subcarriers
NST = 52; %total number of subcarriers
dF = CS/64; %subcarrier frequency spacing
TFFT = 1/dF; % IFFT/FFT period
k = -NST/2:NST/2;

z = zeros(1,3); % length of zero block in short preamble
SPS = 1+1j; % symbol in the short preamble 

SP = sqrt(13/6)*[0 0 SPS z -SPS z SPS z -SPS z -SPS z SPS z 0 z -SPS z -SPS z SPS z SPS z SPS z SPS 0 0];
    
c = sqrt(13/6)*[0 0 1+1j 0 0 0 -1-1j 0 0 0 1+1j 0 0 0 -1-1j 0 0 0 -1-1j 0 0 0 1+1j 0 0 0 0 ...
0 0 0 -1-1j 0 0 0 -1-1j 0 0 0 1+1j 0 0 0 1+1j 0 0 0 1+1j 0 0 0 1+1j 0 0];
% sqrt(13/6) scales so that short preamble has unit average power, get 13/6
% from 52 total subcarriers with only 24 active for the short preamble
%SP = [zeros(1,6) SP zeros(1,6)];

t = (0:63)*(TFFT/64);
myiSP = NaN(length(t),1);
for jj = 1:length(t)
%for ii = 1:length(k)
    myiSP(jj) = sum(SP.*exp(1j*2*pi*k*dF*t(jj)))/64;
    %myiSP(jj) = sum(SP(7:end-6).*exp(1j*2*pi*k*dF*t(jj)));
%end
end
SP2 = [0 SP(28:end) zeros(1,11) SP(1:26)].';
iSP = ifft(SP2,64); % this matches annex L

wt = [0.5;ones(159,1);0.5];
ceSPT = wt.*[myiSP;myiSP;myiSP(1:33)]; %cyclically extended short preamble in time

%%
%tmp = [SP;c].';

LP = [1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 0 ...
1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1];
% the long preamble symbols, already has average unit power

myiLP = NaN(length(t),1);
for jj = 1:length(t)
    myiLP(jj) = sum(LP.*exp(1j*2*pi*k*dF*t(jj)))/64;
end
LP2 = [zeros(1,6) LP zeros(1,5)];
%LP2 = [0 LP(28:end) zeros(1,11) LP(1:26)].';
iLP = ifft(LP2,64).';
