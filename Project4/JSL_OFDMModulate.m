


function [msgTX,msgBin] = JSL_OFDMModulate(rate,numBytes)
% ce = 1; %channel enable
% fD = 100;
% snrdb = 10;
% numBytes = 120;
numBytes2 = numBytes+2;
rateStruct = mcsInfo(rate);
numSym = ceil((16+8*numBytes2+6)/rateStruct.NDBPS); %defined in the Standard
%
msgBin = randi(2,numBytes2,8)-1; 
%+2 bytes for the SERVICE field that would be prepended to DATA

numBits = numel(msgBin)+6+16;% +6 b/c required to have tail of atleast 6 zeros
pad = rateStruct.NDBPS*ceil(numBits/rateStruct.NDBPS)-numBits; % number of pad bits
msgBin = [zeros(16,1); msgBin(:) ; zeros(pad+6,1)]; %+6 for required tail
msgBin = double(msgBin);

encodedSIGNAL = step(rateStruct.hConvEnc,msgBin);
%q = reshape(encodedSIGNAL,rateStruct.NCBPS,numSym);
q = reshape(encodedSIGNAL,rateStruct.NCBPS,[]);
%q = [(1000:1047).' (1048:1095).'];
%
% interleaving the data
s = max([rateStruct.NBPSC/2 1]); %number of coded bits per subcarrier
kint = (0:rateStruct.NCBPS-1).'; % initial index of data in OFDM symbol
iint = (rateStruct.NCBPS/16)*mod(kint,16)+floor(kint/16); % index after 1st interleave
jint = s*floor(iint/s)+mod((iint+rateStruct.NCBPS-floor(16*iint/rateStruct.NCBPS)),s); %index after second interleace

int1 = NaN(rateStruct.NCBPS,numSym);%data after 1st interleave
int2 = NaN(rateStruct.NCBPS,numSym);%data after 2nd interleave

for qq = 1:numSym
    int1(iint+1,qq) = q(kint+1,qq);
    int2(jint+1,qq) = int1(iint+1,qq);
end
sint2 = size(int2);

% modulating the data
uB2 = (reshape(int2,rateStruct.NBPSC,sint2(1)/rateStruct.NBPSC,numSym)); % data reshaped for modulation
uD = NaN(rateStruct.NMSPOS,numSym); % grouped bits converted to decimal
uM = NaN(rateStruct.NMSPOS,numSym); % decimal converted modulation symbol
for nn = 1:numSym
    uD(:,nn) = bi2de(uB2(:,:,nn).');
    uM(:,nn) = step(rateStruct.hMod,uD(:,nn));
end

% OFDM modulate
%numSymp = max([280 numSym]);
LFSR = ones(7,1);
pilotPolar = NaN(2^length(LFSR)-1,1);
for nn = 1:2^length(LFSR)-1
    nxt = mod(LFSR(end)+LFSR(4),2);
    pilotPolar(nn) = nxt;
    LFSR = [nxt; LFSR(1:end-1)];
end
pilotPolar = -2*pilotPolar+1;
pilotPolarInd = mod(0:numSym,127)+1; % used for indexing purposes later
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

v = NaN(53,numSym)+1j*NaN(53,numSym);
v(M+27,:) = uM;
%q(27,:) = DCN;
v(pilotInds,:) = pnn(:,1:numSym);

d = [padL;v;padR];
%
D = ifft(ifftshift(d,1),64);
clear padL padR acp 
D2 = [D(49:64,:);D]; % add cyclic prefix
w = ones(80,1);
w([1 end]) = 0.5; %discrete window function in 802.11
D3 = w.*D2;
a = 80;
u = (0:numSym)*(a-1)+1;
U = [u(1:end-1); u(2:end)];
%
cont = zeros(numSym*(a-1)+1,numSym);
for ii = 1:numSym
    cont(U(1,ii):U(2,ii),ii)=D3(:,ii);
end
msgTX = sum(cont,2);
end