%% ECE 408 - Wireless Communications
% Project 4 - MIMO OFDM
% OFDM part
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
%% setting up mod/demod
% t = [0:1e-3:pi/2 (pi/2)*ones(1,2e3)  pi/2:1e-3:pi];
% w = sin(t).^2;
% plot(1:length(t),w)

rateStruct = mcsInfo(6);

w = ones(80,1);
w([1 end]) = 0.5; %discrete window function in 802.11
%
%tic
MO = 4;
k = log(MO)/log(2);
p2 = 2.^(0:k-1).'; %powers of 2
p = sqrt(MO);
numSym = 5;
binData = randi(2,numSym*48,k)-1;
decData = (binData*p2)+1;
decData = reshape(decData,48,[]);
const = (-(p-1):2:(p-1))+1j*((p-1):-2:-(p-1)).'; 
acp = sum(abs(const).^2,'all')/numel(const);
constV = const(:)/sqrt(acp);%unit power

x = constV(decData);
% msp = mean(abs(rd).^2); %mean symbol power

numSymp = max([280 numSym]);
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
pnn = (pp2.*pilotVals).'; %pilots and nulls

clear numSymp LFSR pilotPolar pilotPolarInd pp2
%
%numSym = 2; %number of OFDM symbols being generated
padL = zeros(6,numSym);
padR = zeros(5,numSym);
%DCN = zeros(1,numSym);

q = NaN(53,numSym)+1j*NaN(53,numSym);
q(M+27,:) = x;
%q(27,:) = DCN;
q(pilotInds,:) = pnn(:,1:numSym);

d = [padL;q;padR];
D = ifft(ifftshift(d,1),64);
clear padL padR acp 
D2 = [D(49:64,:);D]; % add cyclic prefix
D3 = w.*D2;
a = 80;
u = (0:numSym)*(a-1)+1;
U = [u(1:end-1); u(2:end)];
%
cont = zeros(numSym*(a-1)+1,numSym);
for ii = 1:numSym
    cont(U(1,ii):U(2,ii),ii)=D3(:,ii);
end
D4 = sum(cont,2);

h = genRayleighFadingV3(u(end),1,1);

ce = 0; %channel enable
if isequal(ce, 1)
    D5 = h.*D4;
else
    D5 = D4;
end
%
sp = mean(abs(D5).^2);
snrdb = 30;
s2 = sp.*10.^(-snrdb/10);

n = (randn(u(end),1)+1j*randn(u(end),1)).*sqrt(s2/2);
% cv = n'*n/a;
np = randi(10,2,1);
n1 = [];%sqrt(s2/2)*(+randn(np(1),1)+1j*randn(np(1),1));
n2 = [];%sqrt(s2/2)*(+randn(np(2),1)+1j*randn(np(2),1));
r2 = [n1;D5+n;n2];
cte = NaN(length(r2)-77,1); %coarse timing estimate
for ii = 1:length(cte)
    cte(ii) = r2(ii:ii+13)'*r2(ii+77-13:ii+77);
end
%%
r3 = NaN(a,numSym);
for ii = 1:numSym
    r3(:,ii) = r2(u(ii):u(ii+1));
end



%%
r = r2(17:end,:);
R = fftshift(fft(r,64,1),1);

R1 = R(7:59,:);
R1(pilotInds,:) = [];

dm = NaN(48,numSym,MO); %distance metric
for mm = 1:MO
    dm(:,:,mm) = abs(R1-constV(mm));
end
[~,ind] = min(dm,[],3);
binRX = fliplr(dec2bin(ind(:)-1)-'0');
check = sum(abs(binRX-binData),'all');
ber = check/numel(binData)
%toc

%% QAM16
M = 16;
n = sqrt(M);
const = (-n:2:n)+1j*(n:-2:-n).'; 
acp = sum(abs(const).^2,'all')/numel(const);
constV = const(:)/sqrt(acp);%unit power
x = constV(randi(16,48,1)); %randomly generated data
msp = mean(abs(x).^2); %mean symbol power

%% 16PSK
MO = 16;
m = 0:MO-1;
const = exp(1j*2*pi*m/MO);
acp = sum(abs(const).^2)/numel(const);

plot(const,'*')
hold on
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'r')