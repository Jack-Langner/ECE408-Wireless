%% ECE 408 - Wireless Communications
% Project 4 - MIMO OFDM
% Jack Langner - MATLAB 2019b
% Due April 29, 2020

%% System Parameters
dF = 20e6/64;
TFFT = 1/dF;
TG = TFFT/4;
dt = TFFT/64;
t = 0:dt:TFFT-dt;
k = (-26:26).';
k2 = (-32:31).';
r = NaN(64,1)+1j*NaN(64,1);
%% borrowing from myself
numSymp = 280;
LFSR = ones(7,1);
pilotPolar = NaN(2^length(LFSR)-1,1);
for nn = 1:2^length(LFSR)-1
    nxt = mod(LFSR(end)+LFSR(4),2);
    pilotPolar(nn) = nxt;
    LFSR = [nxt; LFSR(1:end-1)];
end
pilotPolar = -2*pilotPolar+1;
pilotPolarInd = mod(0:numSymp-1,127)+1; % used for indexing purposes later
pp2 = pilotPolar(pilotPolarInd);
M = NaN(48,1);%NaN(rateStruct.NMSPOS,1);
LSI = (0:48-1).';%(0:rateStruct.NMSPOS-1).'; %logical subcarrier index
M(LSI>=0 & LSI<=4) = LSI(LSI>=0 & LSI<=4)-26;
M(LSI>=5 & LSI<=17) = LSI(LSI>=5 & LSI<=17)-25;
M(LSI>=18 & LSI<=23) = LSI(LSI>=18 & LSI<=23)-24;
M(LSI>=24 & LSI<=29) = LSI(LSI>=24 & LSI<=29)-23;
M(LSI>=30 & LSI<=42) = LSI(LSI>=30 & LSI<=42)-22;
M(LSI>=43 & LSI<=47) = LSI(LSI>=43 & LSI<=47)-21;
%
%Pk = zeros(53,1);
pilotInds = [6 20 27 34 48]; %includes DC null
pilotVals = [1 1 0 1 -1];
%Pk(pilotInds) = pilotVals;
w = (pp2.*pilotVals).';

clear numSymp LFSR pilotPolar pilotPolarInd pp2
%
numSym = 2; %number of OFDM symbols being generated
padL = zeros(6,numSym);
padR = zeros(5,numSym);
%DCN = zeros(1,numSym);

% QAM16
MO = 16;
n = sqrt(MO)-1;
%n = N-1;
const = (-n:2:n)+1j*(n:-2:-n).'; 
acp = sum(abs(const).^2,'all')/numel(const);
constV = const(:)/sqrt(acp);%unit power
rd = constV(randi(MO,48,numSym)); %randomly generated data
%msp = mean(abs(rd).^2); %mean symbol power

q = NaN(53,numSym)+1j*NaN(53,numSym);
q(M+27,:) = rd;
%q(27,:) = DCN;
q(pilotInds,:) = w(:,1:numSym);

d = [padL;q;padR];
D = ifft(ifftshift(d,1),64);
clear padL padR acp 
D2 = [D(49:64,:);D];

r = D2;%+sqrt(0.0001/2)*(randn(64,numSym)+1j*randn(64,numSym));

%R = fft(fftshift(r));
%%
subplot(2,2,1)
plot(real(D))
title('real(D)')
subplot(2,2,2)
plot(imag(D))
title('imag(D)')
subplot(2,2,3)
plot(real(r))
title('real(r)')
subplot(2,2,4)
plot(imag(r))
title('imag(r)')

%%

Omod = comm.OFDMModulator('FFTLength',64,'NumGuardBandCarriers',[6;5],...
    'InsertDCNull',true,'PilotInputPort',true,'PilotCarrierIndices',[12; 26;40;54],...
    'CyclicPrefixLength',16,'NumSymbols',numSym);

r2 = step(Omod,rd,w([1 2 4 5],1:numSym));

% Demod = comm.OFDMDemodulator(Omod);
% v = step(Demod,r);
% v2 = step(Demod,r2);
%%

subplot(1,2,1)
plot(real([r r2]))
title('real')
legend('r','r_2')
subplot(1,2,2)
plot(imag([r r2]))
title('imag')
legend('r','r_2')
%% generate nakagami rand vars.
pd = makedist('Nakagami','mu',5,'omega',2);
r = random(pd,1e4,1);
mean(r)
histogram(r,100,'Normalization','pdf')

%% QAM16
M = 16;
n = sqrt(M);
const = (-n:2:n)+1j*(n:-2:-n).'; 
acp = sum(abs(const).^2,'all')/numel(const);
constV = const(:)/sqrt(acp);%unit power
rd = constV(randi(16,48,1)); %randomly generated data
msp = mean(abs(rd).^2); %mean symbol power

%% 16PSK
MO = 16;
m = 0:MO-1;
const = exp(1j*2*pi*m/MO);
acp = sum(abs(const).^2)/numel(const);

plot(const,'*')
hold on
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'r')