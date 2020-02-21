
msg = [0x04 0x02 0x00 0x2E 0x00 0x60 0x08 0xCD 0x37 0xA6 ...
 0x00 0x20 0xD6 0x01 0x3C 0xF1 0x00 0x60 0x08 0xAD ...
 0x3B 0xAF 0x00 0x00 0x4A 0x6F 0x79 0x2C 0x20 0x62 ...
 0x72 0x69 0x67 0x68 0x74 0x20 0x73 0x70 0x61 0x72 ...
 0x6B 0x20 0x6F 0x66 0x20 0x64 0x69 0x76 0x69 0x6E ...
 0x69 0x74 0x79 0x2C 0x0A 0x44 0x61 0x75 0x67 0x68 ...
 0x74 0x65 0x72 0x20 0x6F 0x66 0x20 0x45 0x6C 0x79 ...
 0x73 0x69 0x75 0x6D 0x2C 0x0A 0x46 0x69 0x72 0x65 ...
 0x2D 0x69 0x6E 0x73 0x69 0x72 0x65 0x64 0x20 0x77 ...
 0x65 0x20 0x74 0x72 0x65 0x61 0x67 0x33 0x21 0xB6];
%% convolutional encoding the signal frame
% 802.11n-2012 pg 2641/2793, annex L
SIGNAL = [1 0 1 1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0].';

trel = poly2trellis(7 ,[133 171]);

hConvEnc = comm.ConvolutionalEncoder(trel,'TerminationMethod','Truncated');
%encodedSIGNAL = hConvEnc(SIGNAL,1);
encodedSIGNAL1 = step(hConvEnc,SIGNAL);
encodedSIGNAL2 = convenc(SIGNAL,trel);

r = [1 1 0 1 0 0 0 1 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1 1 1 1 1 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0].';
tmp = [encodedSIGNAL1 encodedSIGNAL2 r];

%% interleaving the signal frame
% given in annex L after interleaving
LInter = [1 0 0 1 0 1 0 0 1 1 0 1 0 0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 0 0 1 1 0 0 1 0 0 1 0 0 1 0 0 1 0 1 0 0].';

NCBPS = 48;

s = max([1/2 1]); %number of coded bits per subcarrier
kint = (0:NCBPS-1).';
iint = (NCBPS/16)*mod(kint,16)+floor(kint/16);
jint = s*floor(iint/s)+mod((iint+NCBPS-floor(16*iint/NCBPS)),s);
%[k mink(iint,48) mink(jint,48)]

int1 = NaN(NCBPS,1);
int2 = NaN(NCBPS,1);
%myScrambling;
for nn = 1:NCBPS
    int1(iint(nn)+1) = r(kint(nn)+1);
    int2(jint(nn)+1) = int1(iint(nn)+1);
end

tmp = [LInter int2];

% would scramble data here, but signal frame not scrambled
%% modulating the data
hMod = comm.BPSKModulator('PhaseOffset',pi);
modSIG = step(hMod,int2);

%% OFDM modulate the data
M = NaN(NCBPS,1);

M(kint>=0 & kint<=4) = kint(kint>=0 & kint<=4)-26;
M(kint>=5 & kint<=17) = kint(kint>=5 & kint<=17)-25;
M(kint>=18 & kint<=23) = kint(kint>=18 & kint<=23)-24;
M(kint>=24 & kint<=29) = kint(kint>=24 & kint<=29)-23;
M(kint>=30 & kint<=42) = kint(kint>=30 & kint<=42)-22;
M(kint>=43 & kint<=47) = kint(kint>=43 & kint<=47)-21;
%%

modSigPilot = zeros(NCBPS+5,1);
for nn = 1:NCBPS
    modSigPilot(M(nn)+27) = modSIG(nn);
end
% modSigPilot = [zeros(6,1); modSigPilot; zeros(5,1)];
% modSigPilot(33) = 0;
%%
pilotInd = [12 26 40 54]-6;
pilotVal = [1 1 1 -1];

modSigPilot(pilotInd) = pilotVal;
% 
% txSig = ifft(modSigPilot,64);
% txSig = [txSig; txSig(1:16)];
%%
hOFDMmod = comm.OFDMModulator;%('CyclicPrefixLength',0);%('InsertDCNull',true,'PilotInputPort',true);
txSIG0 = step(hOFDMmod,modSigPilot);

EbNodB = 0:0.5:12;
codeRate = 1/2;
M = 2;
k = log2(M);
SNRdB = EbNodB+10*log10(k*codeRate);
jj = 13;

%noise = sqrt(10^(-SNRdB(jj)/10)/2)*(randn(length(txSIG0),1)+1j*randn(length(txSIG0),1));
%txSIG = txSIG0+noise;
txSIG = awgn(txSIG0,SNRdB(jj),'measured');
hOFDMdemod = comm.OFDMDemodulator;
rxSIG = step(hOFDMdemod,txSIG);
nullCarr = [6 20 27 34 48];
rxSIG(nullCarr) = [];
%% demodulating the data
hDemod = comm.BPSKDemodulator('PhaseOffset',pi);
demodSIG = step(hDemod,rxSIG);
%% deinterleaving the signal frame
ide = s*floor(jint/s)+mod(jint+floor(16*jint/NCBPS),s);
kde = 16*ide-(NCBPS-1)*floor(16*ide/NCBPS);
deint2 = NaN(NCBPS,1);
deint1 = NaN(NCBPS,1);
for nn = 1:NCBPS
    deint2(ide(nn)+1) = demodSIG(jint(nn)+1);
    deint1(kde(nn)+1) = deint2(ide(nn)+1);
end
%% viterbi decoding the signal frame
%hVitDec = comm.ViterbiDecoder(trel,'TracebackDepth',24,'TerminationMethod','Truncated');
hVitDec = comm.ViterbiDecoder(trel,'TracebackDepth',24,'InputFormat','Hard','TerminationMethod','Truncated');
decodedSIGNAL = step(hVitDec, deint1);
decSIG = vitdec(encodedSIGNAL2,trel,24,'term','hard');
tmp2 = [SIGNAL decodedSIGNAL decSIG];