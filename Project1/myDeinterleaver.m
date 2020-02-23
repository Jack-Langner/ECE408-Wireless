
function deint1 = myDeinterleaver(rateStruct,numSym,binMat)
% deinterleaving the signal frame
s = max([rateStruct.NBPSC/2 1]); %number of coded bits per subcarrier
kint = (0:rateStruct.NCBPS-1).'; % initial index of data in OFDM symbol
iint = (rateStruct.NCBPS/16)*mod(kint,16)+floor(kint/16); % index after 1st interleave
jint = s*floor(iint/s)+mod((iint+rateStruct.NCBPS-floor(16*iint/rateStruct.NCBPS)),s); %index after second interleave
ide = s*floor(jint/s)+mod(jint+floor(16*jint/rateStruct.NCBPS),s);
kde = 16*ide-(rateStruct.NCBPS-1)*floor(16*ide/rateStruct.NCBPS);
% ide and kde are the indices after the 1st and 2nd deinterleaving,
% respectively
deint2 = NaN(rateStruct.NCBPS,numSym); % data after 1st deinterleave
deint1 = NaN(rateStruct.NCBPS,numSym); % data after 2nd deinterleave
for qq = 1:numSym
for nn = 1:rateStruct.NCBPS
    deint2(ide(nn)+1,qq) = binMat(jint(nn)+1,qq);
    deint1(kde(nn)+1,qq) = deint2(ide(nn)+1,qq);
end
end

end