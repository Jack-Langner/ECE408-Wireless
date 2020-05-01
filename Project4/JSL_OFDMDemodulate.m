

function decodedSIGNAL = JSL_OFDMDemodulate(rate,numBytes,msgRX)
numBytes2 = numBytes+2;
rateStruct = mcsInfo(rate);
numSym = ceil((16+8*numBytes2+6)/rateStruct.NDBPS); %defined in the Standard
a = 80;
u = (0:numSym)*(a-1)+1;

r3 = NaN(a,numSym)+1j*NaN(a,numSym);
for qq = 1:numSym
    r3(:,qq) = msgRX(u(qq):u(qq+1));
end
r3([1 80],:) = r3([65 16],:); % replacing the overlap samples

r = r3(17:end,:);
R = fftshift(fft(r,64,1),1);

R1 = R(7:59,:);
pilotInds = [6 20 27 34 48]; %includes DC null
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
s = max([rateStruct.NBPSC/2 1]); %number of coded bits per subcarrier
kint = (0:rateStruct.NCBPS-1).'; % initial index of data in OFDM symbol
iint = (rateStruct.NCBPS/16)*mod(kint,16)+floor(kint/16); % index after 1st interleave
jint = s*floor(iint/s)+mod((iint+rateStruct.NCBPS-floor(16*iint/rateStruct.NCBPS)),s); %index after second interleace
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
end