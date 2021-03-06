
% Right off the bat, I made a conscious decision to not make this a
% function becasue that would require more effort I feel I completed the
% project as intended at this point and would like to move on. 
% right now the script is set to read in text from a file named tmp.txt in
% the current working directory, feel free to change the contents of that
% file or the file that is read from. The SNR is controlled with the
% variable jj which is defined on line 43. yy controls the data rate, low
% rates give reliable TX/RX, high rates not so much. I also recognize that
% I may have started to skimp on comments and came up with odd variable
% names, sorry it happens after a while. 
% Have fun. JSL
clear;clc

modes = [6 9 12 18 24 36 48 54];
yy = 2;
txRS = mcsInfo(modes(yy));

crcPoly= [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
genCRC= comm.CRCGenerator(crcPoly);
%
fileID = fopen('tmp.txt','r');
[msg, numBytes] = fscanf(fileID,'%c');
msgASCII = (double(msg)).';
msgBin = de2bi(msgASCII,8);

msgBin = ([randi(2,10,8)-1;msgBin]).'; % adding minimum MAC header length
msgBin = msgBin(:);
msgBin = genCRC(msgBin); %this is used to determine numSym
numSym = ceil((16+length(msgBin)+6)/txRS.NDBPS); %defined in the Standard
cBin = de2bi(length(msgBin)/8,12);

numBits = numel(msgBin)+6+16;% +6 b/c required to have tail of atleast 6 zeros
pad = txRS.NDBPS*ceil(numBits/txRS.NDBPS)-numBits; % number of pad bits
msgBin = [zeros(16,1); msgBin(:); zeros(pad+6,1)]; %+6 for required tail
msgBin = double(msgBin);
%%
EbNodB = (0:0.5:12).';
lenEbNo = length(EbNodB);
%codeRate = 1/2;
%M = 2;
%k = log2(M);
%jj = 15;
BERsim = zeros(lenEbNo,length(modes));
SNRdB = EbNodB+10*log10(txRS.NBPSC*txRS.eccRate);
lenSNR = length(SNRdB);
jj = 8;
% wb = waitbar(0,'Beginning simulation');
% tstart = clock;
% for yy = 1:length(modes)
sigRS = mcsInfo(6); % rate struct for the signal frame
sigRS.hVitDec.TracebackDepth = 18; % have to change for SIGframe b/c small msg
SIGframe = [txRS.R1R4 0 cBin mod(sum([txRS.R1R4 cBin]),2) zeros(1,6)].';
encSIGframe = step(sigRS.hConvEnc, SIGframe);


intLSIGframe= myInterleaver(sigRS,1,encSIGframe);
modSIGframe = step(sigRS.hMod,intLSIGframe);

LFSR = ones(7,1);
pilotPolar = NaN(2^length(LFSR)-1,1);
for nn = 1:2^length(LFSR)-1
    nxt = mod(LFSR(end)+LFSR(4),2);
    pilotPolar(nn) = nxt;
    LFSR = [nxt; LFSR(1:end-1)];
end
clear LFSR
pilotPolar = -2*pilotPolar+1;
pilotPolarInd = mod(0:numSym,127)+1; % used for indexing purposes later

%%
initState = randi(2,7,1)-1; % random bit sequence
%initializer the scrambler in a random state per 802.11
LFSR = initState; 
seq = NaN(length(msgBin),1); 
% sequence generated by the LSFR for scrambling
for nn = 1:length(seq)
    nxt = mod(LFSR(end)+LFSR(4),2); % newest bit into LFSR
    seq(nn) = nxt;
    LFSR = [nxt; LFSR(1:end-1)];
end

msgScr = mod(double(msgBin)+seq,2); 
% the scrambled message, xor of data and scrambler sequence

% convolutional encoding the data frame for the given data rate
encodedSIGNAL = step(txRS.hConvEnc,msgScr);
encSIGMat = reshape(encodedSIGNAL,txRS.NCBPS,numSym);
% reshape encSig b/c everything that follows is based on size of OFDM frame

int2 = myInterleaver(txRS,numSym,encSIGMat);
sint2 = size(int2);

% modulating the data
uB2 = (reshape(int2,txRS.NBPSC,sint2(1)/txRS.NBPSC,numSym)); % data reshaped for modulation
uD = NaN(txRS.NMSPOS,numSym); % grouped bits converted to decimal
uM = NaN(txRS.NMSPOS,numSym); % decimal converted modulation symbol
for nn = 1:numSym
    uD(:,nn) = bi2de(uB2(:,:,nn).');
    uM(:,nn) = step(txRS.hMod,uD(:,nn));
end

PL = [modSIGframe uM]; % complete payload ti be sent
%%
% M is a function in 802.11 that maps LSI to OFDM tones
M = NaN(sigRS.NMSPOS,1);
LSI = (0:sigRS.NMSPOS-1).'; %logical subcarrier index
M(LSI>=0 & LSI<=4) = LSI(LSI>=0 & LSI<=4)-26;
M(LSI>=5 & LSI<=17) = LSI(LSI>=5 & LSI<=17)-25;
M(LSI>=18 & LSI<=23) = LSI(LSI>=18 & LSI<=23)-24;
M(LSI>=24 & LSI<=29) = LSI(LSI>=24 & LSI<=29)-23;
M(LSI>=30 & LSI<=42) = LSI(LSI>=30 & LSI<=42)-22;
M(LSI>=43 & LSI<=47) = LSI(LSI>=43 & LSI<=47)-21;
%
pilotInd = [12 26 40 54]-6; % indeces where pilot tones are inserted
pilotVal = [1 1 1 -1]; % symbols on the pilot tones
modSigPilot = zeros(sigRS.NMSPOS+5,numSym+1); 
% modulated signal with pilot tones inserted and DC null
for qq = 1:numSym+1
    for nn = 1:sigRS.NMSPOS
        modSigPilot(M(nn)+27,qq) = PL(nn,qq); %indexing game
    end
    modSigPilot(pilotInd,qq) = pilotVal*pilotPolar(pilotPolarInd(qq));
end

% OFDM modulate the data
% using vanilla OFDM modulator from MathWorks because I couldnt figure out
% how 802.11 implements IFFT
txSIG = NaN(sigRS.hOFDMmod.FFTLength+sigRS.hOFDMmod.CyclicPrefixLength,numSym+1);
for qq = 1:numSym+1
txSIG(:,qq) = step(sigRS.hOFDMmod,modSigPilot(:,qq));
end
%%
txSIG = awgn(txSIG,SNRdB(jj),'measured'); %adding noise

% can use the sigRS struct becasue the OFDM portion of TR/RX is the same
% regardless of MCS chosen, so once the SIGNAL frame is processed the rest
% of the message can be processed with only knowledge gained from the
% SIGNAL frame and assuming the size of the MAC header
rxSIG = NaN(sigRS.NMSPOS,numSym+1); % received OFDM symbols
nullCarr = [6 20 27 34 48]; % no data symbols at these indices
for qq = 1:numSym+1
tmp = step(sigRS.hOFDMdemod,txSIG(:,qq));
tmp(nullCarr) = [];
rxSIG(:,qq) = tmp;
end

rxSIGframe = rxSIG(:,1);
rxSIG = rxSIG(:,2:end);

demodSIGframe = step(sigRS.hDemod,rxSIGframe);
deintSIGframe = myDeinterleaver(sigRS,1,demodSIGframe);
decodSIGframe = step(sigRS.hVitDec,deintSIGframe);
decodSIGframe = decodSIGframe.';
%%
rxRateBin = decodSIGframe(1:4);
rxLenBin = decodSIGframe(6:17);
rxRate = miniMap(rxRateBin);
rxLen = bi2de(rxLenBin);

rxRS = mcsInfo(rxRate);
numRxSym = ceil((16+8*rxLen+6)/rxRS.NDBPS); %defined in the Standard
totBitsRx = numRxSym*48*rxRS.NBPSC*rxRS.eccRate;
numPadBits = totBitsRx-8*rxLen-6-16;
numBitRx = totBitsRx-16-6-numPadBits-8*14; %assuming min MAC header
% data rate is 36 Mbps
%%
% demodulating the data
demodSIG = NaN(rxRS.NMSPOS,numRxSym);
uB3 = NaN(rxRS.NBPSC,rxRS.NMSPOS,numRxSym);
uB4 = NaN(rxRS.NMSPOS*rxRS.NBPSC,numRxSym);
% just indexing/transpose games to get demodulated symbols
% into correct bit pattern
for qq = 1:numRxSym
    demodSIG(:,qq) = step(rxRS.hDemod,rxSIG(:,qq));
    uB3(:,:,qq) = (de2bi(demodSIG(:,qq),rxRS.NBPSC)).';
    tmp = uB3(:,:,qq);
    uB4(:,qq) = tmp(:);
end

deintRxMsg = myDeinterleaver(rxRS,numRxSym,uB4);
rxencMsg = deintRxMsg(:); % received msg, still encoded
% viterbi decoding the signal frame
decodedSIGNAL = step(rxRS.hVitDec, rxencMsg);
rxmsgBin = mod(decodedSIGNAL+seq,2);
% descramble the decoded msg, xor with same pattern recovers original 

BER = sum(abs(rxmsgBin-msgBin))/length(msgBin);

QQQ = rxmsgBin(16+1:end-(numPadBits+6)); % 'processing' the SERVICE field and pad bits
QQQrs = (reshape(QQQ,8,length(QQQ)/8)).';% recover matrix structure for ASCII conversion
dataBin = QQQrs(11:end-4,:); % 'processing' the MAC header and CRC32
% Note, the CRC32 would be used as final check to see if any bit errors had
% occured at this point
dataDec = bi2de(dataBin);
promisedLand = convertCharsToStrings(char(dataDec))