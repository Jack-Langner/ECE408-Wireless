
function BER = berSim80211(numBytes,rateStruct,numSym,SNRdB)
LFSR = ones(7,1);
pilotPolar = NaN(2^length(LFSR)-1,1);
for nn = 1:2^length(LFSR)-1
    nxt = mod(LFSR(end)+LFSR(4),2);
    pilotPolar(nn) = nxt;
    LFSR = [nxt; LFSR(1:end-1)];
end
clear LFSR
pilotPolar = -2*pilotPolar+1;
pilotPolarInd = mod(0:numSym-1,127)+1;

msgBin = randi(2,numBytes+2,8)-1;
%
numBits = numel(msgBin)+6;% +6 b/c required to have tail of atleast 6 zeros
pad = rateStruct.NDBPS*ceil(numBits/rateStruct.NDBPS)-numBits; % number of pad bits
msgBin = [msgBin(:) ; zeros(pad+6,1)]; %+6 for required tail
msgBin = double(msgBin);
%seq = myScrambling([1 0 1 1 1 0 1].',length(msgBin));
initState = randi(2,7,1)-1;
LFSR = initState;%[1 0 1 1 1 0 1].';%ones(7,1); 
seq = NaN(length(msgBin),1);

for nn = 1:length(seq)
    nxt = mod(LFSR(end)+LFSR(4),2);
    seq(nn) = nxt;
    LFSR = [nxt; LFSR(1:end-1)];
end

msgScr = mod(double(msgBin)+seq,2);
%% convolutional encoding the data frame

encodedSIGNAL = step(rateStruct.hConvEnc,msgScr);
q = reshape(encodedSIGNAL,rateStruct.NCBPS,numSym);
%% interleaving the data

s = max([rateStruct.NBPSC/2 1]); %number of coded bits per subcarrier
kint = (0:rateStruct.NCBPS-1).';
iint = (rateStruct.NCBPS/16)*mod(kint,16)+floor(kint/16);
jint = s*floor(iint/s)+mod((iint+rateStruct.NCBPS-floor(16*iint/rateStruct.NCBPS)),s);
%%
int1 = NaN(rateStruct.NCBPS,numSym);
int2 = NaN(rateStruct.NCBPS,numSym);

for qq = 1:numSym
for nn = 1:rateStruct.NCBPS
    int1(iint(nn)+1,qq) = q(kint(nn)+1,qq);
    int2(jint(nn)+1,qq) = int1(iint(nn)+1,qq);
end
end
sint2 = size(int2);

%% modulating the data

uB2 = (reshape(int2,rateStruct.NBPSC,sint2(1)/rateStruct.NBPSC,numSym)); % data reshaped for modulation
uD = NaN(rateStruct.NMSPOS,numSym);
uM = NaN(rateStruct.NMSPOS,numSym);
for nn = 1:numSym
    uD(:,nn) = bi2de(uB2(:,:,nn).');
    uM(:,nn) = step(rateStruct.hMod,uD(:,nn));
end

%% OFDM modulate the data
M = NaN(rateStruct.NMSPOS,1);
LSI = (0:rateStruct.NMSPOS-1).'; %logical subcarrier index
M(LSI>=0 & LSI<=4) = LSI(LSI>=0 & LSI<=4)-26;
M(LSI>=5 & LSI<=17) = LSI(LSI>=5 & LSI<=17)-25;
M(LSI>=18 & LSI<=23) = LSI(LSI>=18 & LSI<=23)-24;
M(LSI>=24 & LSI<=29) = LSI(LSI>=24 & LSI<=29)-23;
M(LSI>=30 & LSI<=42) = LSI(LSI>=30 & LSI<=42)-22;
M(LSI>=43 & LSI<=47) = LSI(LSI>=43 & LSI<=47)-21;
%%
pilotInd = [12 26 40 54]-6;
pilotVal = [1 1 1 -1];
modSigPilot = zeros(rateStruct.NMSPOS+5,numSym);
for qq = 1:numSym
    for nn = 1:rateStruct.NMSPOS
        modSigPilot(M(nn)+27,qq) = uM(nn,qq);
    end
    modSigPilot(pilotInd,qq) = pilotVal*pilotPolar(pilotPolarInd(qq));
end

%% 
%hOFDMmod = comm.OFDMModulator;%('CyclicPrefixLength',0);%('InsertDCNull',true,'PilotInputPort',true);
txSIG = NaN(rateStruct.hOFDMmod.FFTLength+rateStruct.hOFDMmod.CyclicPrefixLength,numSym);
for qq = 1:numSym
txSIG(:,qq) = step(rateStruct.hOFDMmod,modSigPilot(:,qq));
end

%%

txSIG = awgn(txSIG,SNRdB,'measured');
%hOFDMdemod = comm.OFDMDemodulator;
rxSIG = NaN(rateStruct.NMSPOS,numSym);
nullCarr = [6 20 27 34 48];
for qq = 1:numSym
tmp = step(rateStruct.hOFDMdemod,txSIG(:,qq));
tmp(nullCarr) = [];
rxSIG(:,qq) = tmp;
end
%rxSIG(nullCarr) = [];
%% demodulating the data
demodSIG = NaN(rateStruct.NMSPOS,numSym);
uB3 = NaN(rateStruct.NBPSC,rateStruct.NMSPOS,numSym);
uB4 = NaN(rateStruct.NMSPOS*rateStruct.NBPSC,numSym);
for qq = 1:numSym
    demodSIG(:,qq) = step(rateStruct.hDemod,rxSIG(:,qq));
    uB3(:,:,qq) = (de2bi(demodSIG(:,qq),rateStruct.NBPSC)).';
    tmp = uB3(:,:,qq);
    uB4(:,qq) = tmp(:);
end

%% deinterleaving the signal frame
ide = s*floor(jint/s)+mod(jint+floor(16*jint/rateStruct.NCBPS),s);
kde = 16*ide-(rateStruct.NCBPS-1)*floor(16*ide/rateStruct.NCBPS);

deint2 = NaN(rateStruct.NCBPS,numSym);
deint1 = NaN(rateStruct.NCBPS,numSym);
for qq = 1:numSym
for nn = 1:rateStruct.NCBPS
    deint2(ide(nn)+1,qq) = uB4(jint(nn)+1,qq);
    deint1(kde(nn)+1,qq) = deint2(ide(nn)+1,qq);
end
end
rxencMsg = deint1(:);
%% viterbi decoding the signal frame
decodedSIGNAL = step(rateStruct.hVitDec, rxencMsg);
rxmsgBin = mod(decodedSIGNAL+seq,2);

BER = sum(abs(rxmsgBin-msgBin))/length(msgBin);
end