

function BER = OFDMfading(rate,numBytes,snrdb,ce,fD)
%ce = 1; %channel enable
%fD = 100;
%snrdb = 10;
%numBytes = 120;
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
D4 = sum(cont,2);

h = genRayleighFadingV3(u(end),fD,1);

if isequal(ce, 1)
    D5 = h.*D4;
else
    D5 = D4;
end
%
sp = mean(abs(D5).^2);
s2 = sp.*10.^(-snrdb/10);

n = (randn(u(end),1)+1j*randn(u(end),1)).*sqrt(s2/2);
% cv = n'*n/a;
np = randi(10,2,1);
n1 = [];%sqrt(s2/2)*(+randn(np(1),1)+1j*randn(np(1),1));
n2 = [];%sqrt(s2/2)*(+randn(np(2),1)+1j*randn(np(2),1));
r2 = [n1;D5+n;n2]; % timing offset
% perform timing offset
r3 = NaN(a,numSym)+1j*NaN(a,numSym);
for qq = 1:numSym
    r3(:,qq) = r2(u(qq):u(qq+1));
end
r3([1 80],:) = r3([65 16],:); % replacing the overlap samples

r = r3(17:end,:);
R = fftshift(fft(r,64,1),1);

R1 = R(7:59,:);
R1(pilotInds,:) = [];

% demodulating the data
demodSIG = NaN(rateStruct.NMSPOS,numSym);
uB3 = NaN(rateStruct.NBPSC,rateStruct.NMSPOS,numSym);
uB4 = NaN(rateStruct.NMSPOS*rateStruct.NBPSC,numSym);
% just indexing/transpose games to get demodulated symbols
% into correct bit pattern
for qq = 1:numSym
    demodSIG(:,qq) = step(rateStruct.hDemod,R1(:,qq));
    uB3(:,:,qq) = (de2bi(demodSIG(:,qq),rateStruct.NBPSC)).';
    tmp = uB3(:,:,qq);
    uB4(:,qq) = tmp(:);
end

% deinterleaving the signal frame
ide = s*floor(jint/s)+mod(jint+floor(16*jint/rateStruct.NCBPS),s);
kde = 16*ide-(rateStruct.NCBPS-1)*floor(16*ide/rateStruct.NCBPS);
% ide and kde are the indices after the 1st and 2nd deinterleaving,
% respectively
deint2 = NaN(rateStruct.NCBPS,numSym); % data after 1st deinterleave
deint1 = NaN(rateStruct.NCBPS,numSym); % data after 2nd deinterleave
for qq = 1:numSym
%for nn = 1:rateStruct.NCBPS
    deint2(ide+1,qq) = uB4(jint+1,qq);
    deint1(kde+1,qq) = deint2(ide+1,qq);
%end
end

rxencMsg = deint1(:); % received msg, still encoded
% viterbi decoding the signal frame
decodedSIGNAL = step(rateStruct.hVitDec, rxencMsg);

BER = sum(abs(decodedSIGNAL-msgBin))/length(msgBin);
end