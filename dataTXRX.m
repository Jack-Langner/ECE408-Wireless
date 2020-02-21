% data rate is 36 Mbps
modes = [13 15 5 7 9 11 1 3];

numBytes = 1024;


EbNodB = (0:2:12).';
lenEbNo = length(EbNodB);
%codeRate = 1/2;
%M = 2;
%k = log2(M);
%jj = 15;
BERsim = zeros(lenEbNo,length(modes));
LFSR = ones(7,1);
pilotPolar = NaN(2^length(LFSR)-1,1);
for nn = 1:2^length(LFSR)-1
    nxt = mod(LFSR(end)+LFSR(4),2);
    pilotPolar(nn) = nxt;
    LFSR = [nxt; LFSR(1:end-1)];
end
clear LFSR
pilotPolar = -2*pilotPolar+1;

numIter = 10;
wb = waitbar(0,'Performing simulation');
tstart = clock;
for yy = 1:length(modes)
RS = mcsInfo(modes(yy));
numSym = ceil((16+8*numBytes+6)/RS.NDBPS); %defined in the Standarad
SNRdB = EbNodB+10*log10(RS.NBPSC*RS.eccRate);
pilotPolarInd = mod(0:numSym-1,127)+1;
for ii = 1:numIter
for jj = 1:lenEbNo
    waitbar(((yy-1)*numIter*lenEbNo+(ii-1)*lenEbNo+jj)/(lenEbNo*numIter*length(modes)),wb)
%%
% msg = [0x04 0x02 0x00 0x2E 0x00 0x60 0x08 0xCD 0x37 0xA6 ...
%  0x00 0x20 0xD6 0x01 0x3C 0xF1 0x00 0x60 0x08 0xAD ...
%  0x3B 0xAF 0x00 0x00 0x4A 0x6F 0x79 0x2C 0x20 0x62 ...
%  0x72 0x69 0x67 0x68 0x74 0x20 0x73 0x70 0x61 0x72 ...
%  0x6B 0x20 0x6F 0x66 0x20 0x64 0x69 0x76 0x69 0x6E ...
%  0x69 0x74 0x79 0x2C 0x0A 0x44 0x61 0x75 0x67 0x68 ...
%  0x74 0x65 0x72 0x20 0x6F 0x66 0x20 0x45 0x6C 0x79 ...
%  0x73 0x69 0x75 0x6D 0x2C 0x0A 0x46 0x69 0x72 0x65 ...
%  0x2D 0x69 0x6E 0x73 0x69 0x72 0x65 0x64 0x20 0x77 ...
%  0x65 0x20 0x74 0x72 0x65 0x61 0x67 0x33 0x21 0xB6];
%msgBin = [zeros(2,8); de2bi(msg)].'; %prepend the SERVICE field

% first two rows are the SERVICE field
msgBin = randi(2,numBytes+2,8)-1;
%
numBits = numel(msgBin)+6;% +6 b/c required to have tail of atleast 6 zeros
pad = RS.NDBPS*ceil(numBits/RS.NDBPS)-numBits; % number of pad bits
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

encodedSIGNAL = step(RS.hConvEnc,msgScr);
q = reshape(encodedSIGNAL,RS.NCBPS,numSym);
%% interleaving the data

s = max([RS.NBPSC/2 1]); %number of coded bits per subcarrier
kint = (0:RS.NCBPS-1).';
iint = (RS.NCBPS/16)*mod(kint,16)+floor(kint/16);
jint = s*floor(iint/s)+mod((iint+RS.NCBPS-floor(16*iint/RS.NCBPS)),s);
%tt = [kint iint jint];
%%
int1 = NaN(RS.NCBPS,numSym);
int2 = NaN(RS.NCBPS,numSym);

for qq = 1:numSym
for nn = 1:RS.NCBPS
    int1(iint(nn)+1,qq) = q(kint(nn)+1,qq);
    int2(jint(nn)+1,qq) = int1(iint(nn)+1,qq);
end
end
sint2 = size(int2);

%% modulating the data

uB2 = (reshape(int2,RS.NBPSC,sint2(1)/RS.NBPSC,numSym)); % data reshaped for modulation
uD = NaN(RS.NMSPOS,numSym);
uM = NaN(RS.NMSPOS,numSym);
for nn = 1:numSym
    uD(:,nn) = bi2de(uB2(:,:,nn).');
    uM(:,nn) = step(RS.hMod,uD(:,nn));
end

%% OFDM modulate the data
M = NaN(RS.NMSPOS,1);
LSI = (0:RS.NMSPOS-1).'; %logical subcarrier index
M(LSI>=0 & LSI<=4) = LSI(LSI>=0 & LSI<=4)-26;
M(LSI>=5 & LSI<=17) = LSI(LSI>=5 & LSI<=17)-25;
M(LSI>=18 & LSI<=23) = LSI(LSI>=18 & LSI<=23)-24;
M(LSI>=24 & LSI<=29) = LSI(LSI>=24 & LSI<=29)-23;
M(LSI>=30 & LSI<=42) = LSI(LSI>=30 & LSI<=42)-22;
M(LSI>=43 & LSI<=47) = LSI(LSI>=43 & LSI<=47)-21;
%%
pilotInd = [12 26 40 54]-6;
pilotVal = [1 1 1 -1];
modSigPilot = zeros(RS.NMSPOS+5,numSym);
for qq = 1:numSym
    for nn = 1:RS.NMSPOS
        modSigPilot(M(nn)+27,qq) = uM(nn,qq);
    end
    modSigPilot(pilotInd,qq) = pilotVal*pilotPolar(pilotPolarInd(qq));
end

%% 
%hOFDMmod = comm.OFDMModulator;%('CyclicPrefixLength',0);%('InsertDCNull',true,'PilotInputPort',true);
txSIG = NaN(RS.hOFDMmod.FFTLength+RS.hOFDMmod.CyclicPrefixLength,numSym);
for qq = 1:numSym
txSIG(:,qq) = step(RS.hOFDMmod,modSigPilot(:,qq));
end

%%
%noise = sqrt(10^(-SNRdB(jj)/10)/2)*(randn(length(txSIG0),1)+1j*randn(length(txSIG0),1));
%txSIG = txSIG0+noise;
txSIG = awgn(txSIG,SNRdB(jj),'measured');
%hOFDMdemod = comm.OFDMDemodulator;
rxSIG = NaN(RS.NMSPOS,numSym);
nullCarr = [6 20 27 34 48];
for qq = 1:numSym
tmp = step(RS.hOFDMdemod,txSIG(:,qq));
tmp(nullCarr) = [];
rxSIG(:,qq) = tmp;
end
%rxSIG(nullCarr) = [];
%% demodulating the data
demodSIG = NaN(RS.NMSPOS,numSym);
uB3 = NaN(RS.NBPSC,RS.NMSPOS,numSym);
uB4 = NaN(RS.NMSPOS*RS.NBPSC,numSym);
for qq = 1:numSym
    demodSIG(:,qq) = step(RS.hDemod,rxSIG(:,qq));
    uB3(:,:,qq) = (de2bi(demodSIG(:,qq))).';
    tmp = uB3(:,:,qq);
    uB4(:,qq) = tmp(:);
end

%% deinterleaving the signal frame
ide = s*floor(jint/s)+mod(jint+floor(16*jint/RS.NCBPS),s);
kde = 16*ide-(RS.NCBPS-1)*floor(16*ide/RS.NCBPS);

deint2 = NaN(RS.NCBPS,numSym);
deint1 = NaN(RS.NCBPS,numSym);
for qq = 1:numSym
for nn = 1:RS.NCBPS
    deint2(ide(nn)+1,qq) = uB4(jint(nn)+1,qq);
    deint1(kde(nn)+1,qq) = deint2(ide(nn)+1,qq);
end
end
rxencMsg = deint1(:);
%% viterbi decoding the signal frame
decodedSIGNAL = step(RS.hVitDec, rxencMsg);
rxmsgBin = mod(decodedSIGNAL+seq,2);

BERsim(jj,yy) = BERsim(jj,yy) + sum(abs(rxmsgBin-msgBin))/length(msgBin);
end
end

end
tend = clock;
close(wb)
fprintf('Time for simulation was %.2f seconds.\n',etime(tend,tstart))
BERsim = BERsim/numIter;

semilogy(EbNodB,BERsim)
ttl = ['Bit Error Rate for 802.11-2012,' num2str(numIter) ' Iterations. ' num2str(etime(tend,tstart)) ' seconds to sim.'];
title(ttl)
grid on