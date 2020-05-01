%% ECE 408 - Wireless Communications
% Project 4 - MIMO OFDM
% OFDM part
% Jack Langner - MATLAB 2019b
% Due April 29, 2020

%% setting up mod/demod

numIter = 20;
snrdb = 10:20;

MO = 4;
k = sqrt(MO);
n = 2.^(0:k-1).';
numSym = 5000;
const = (-(k-1):2:(k-1))+1j*((k-1):-2:-(k-1)).'; 
acp = sum(abs(const).^2,'all')/numel(const);
constV = const(:)/sqrt(acp);%unit power

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
pilotInds = [6 20 27 34 48]; %includes DC null
pilotVals = [1 1 0 1 -1];
w = (pp2.*pilotVals).';

clear numSymp LFSR pilotPolar pilotPolarInd pp2
padL = zeros(6,numSym);
padR = zeros(5,numSym);
a = 80;
BER = NaN(numIter,length(snrdb));
wb = waitbar(0,'Beginning simulation');
tstart = clock;
for ii = 1:numIter
    for jj = 1:length(snrdb)
        waitbar((numIter*(ii-1)+jj)/(numIter*length(snrdb)),wb,'Performing Sim');

binData = randi(2,numSym*48,k)-1;
decData = (binData*n)+1;
decData = reshape(decData,48,[]);

rd = constV(decData);
% msp = mean(abs(rd).^2); %mean symbol power


q = NaN(53,numSym)+1j*NaN(53,numSym);
q(M+27,:) = rd;
%q(27,:) = DCN;
q(pilotInds,:) = w(:,1:numSym);

d = [padL;q;padR];
D = ifft(ifftshift(d,1),64);
%clear padL padR acp 
D2 = [D(49:64,:);D];
%
sp = mean(abs(D2).^2);
%snrdb = 20;
s2 = sp.*10.^(-snrdb(jj)/10);

noise = (randn(a,numSym)+1j*randn(a,numSym)).*sqrt(s2/2);
%cv = noise'*noise/a;

r2 = D2+noise;
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
BER(ii,jj) = sum(abs(binRX-binData),'all')/numel(binData);
    end
end
tend = clock;
close(wb)
fprintf('Time for simulation was %.2f seconds.\n',etime(tend,tstart))
BER = mean(BER);
%
BERbpsk = berawgn(snrdb,'psk',2,'nondiff');
semilogy(snrdb,BER,'b-*')
hold on
semilogy(snrdb,BERbpsk,'r-s')
legend('OFDM','bpsk')