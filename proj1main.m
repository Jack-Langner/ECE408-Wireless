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

numIter = 5;
wb = waitbar(0,'Beginning simulation');
tstart = clock;
for yy = 1:length(modes)
RS = mcsInfo(modes(yy));
numSym = ceil((16+8*numBytes+6)/RS.NDBPS); %defined in the Standarad
SNRdB = EbNodB+10*log10(RS.NBPSC*RS.eccRate);
%pilotPolarInd = mod(0:numSym-1,127)+1;
for ii = 1:numIter
for jj = 1:lenEbNo
    newMsg = ['Simulating ' RS.ttl ' data rate'];
    waitbar(((yy-1)*numIter*lenEbNo+(ii-1)*lenEbNo+jj)/(lenEbNo*numIter*length(modes)),wb,newMsg)
    BERsim(jj,yy) = BERsim(jj,yy)+berSim80211(numBytes,RS,numSym,SNRdB(jj));
end
end
end

tend = clock;
close(wb)
fprintf('Time for simulation was %.2f seconds.\n',etime(tend,tstart))
BERsim = BERsim/numIter;
figure
semilogy(EbNodB,BERsim)
ttl = ['Bit Error Rate for 802.11-2012,' num2str(numIter) ' Iterations. ' num2str(etime(tend,tstart)) ' seconds to sim.'];
title(ttl)
grid on