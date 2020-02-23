
function intlv = myInterleaver(rateStruct, numSym, encSIGMat)
s = max([rateStruct.NBPSC/2 1]); %number of coded bits per subcarrier
kint = (0:rateStruct.NCBPS-1).'; % initial index of data in OFDM symbol
iint = (rateStruct.NCBPS/16)*mod(kint,16)+floor(kint/16); % index after 1st interleave
jint = s*floor(iint/s)+mod((iint+rateStruct.NCBPS-floor(16*iint/rateStruct.NCBPS)),s); %index after second interleave
%
int1 = NaN(rateStruct.NCBPS,numSym);%data after 1st interleave
intlv = NaN(rateStruct.NCBPS,numSym);%data after 2nd interleave

for qq = 1:numSym
for nn = 1:rateStruct.NCBPS
    int1(iint(nn)+1,qq) = encSIGMat(kint(nn)+1,qq);
    intlv(jint(nn)+1,qq) = int1(iint(nn)+1,qq);
end
end
end