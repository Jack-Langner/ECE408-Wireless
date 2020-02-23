% data rate is 36 Mbps
modes = [6 9 12 18 24 36 48 54];
numBytes = 1*1024;
EbNodB = (0:0.5:12).';
lenEbNo = length(EbNodB);
%codeRate = 1/2;
%M = 2;
%k = log2(M);
%jj = 15;
BERsim = zeros(lenEbNo,length(modes));
SNRdB = 0:12;
lenSNR = length(SNRdB);
numIter = 10;
wb = waitbar(0,'Beginning simulation');
tstart = clock;
for yy = 1:length(modes)
RS = mcsInfo(modes(yy));
numSym = ceil((16+8*numBytes+6)/RS.NDBPS); %defined in the Standard
%SNRdB = EbNodB+10*log10(RS.NBPSC*RS.eccRate);
%pilotPolarInd = mod(0:numSym-1,127)+1;
for ii = 1:numIter
for jj = 1:lenSNR
    newMsg = ['Simulating ' RS.ttl ' data rate at SNR of ' num2str(SNRdB(jj)) 'dB, iteration # ' num2str(ii)];
    waitbar(((yy-1)*numIter*lenEbNo+(ii-1)*lenEbNo+jj)/(lenEbNo*numIter*length(modes)),wb,newMsg)
    BERsim(jj,yy) = BERsim(jj,yy)+berSim80211(numBytes,RS,numSym,SNRdB(jj));
end
end
end

tend = clock;
close(wb)
fprintf('Time for simulation was %.2f seconds.\n',etime(tend,tstart))
BERsim = BERsim/numIter;
%
figure
semilogy(EbNodB,BERsim,'LineWidth',3)
ttl = ['Bit Error Rate for 802.11-2012 with ' num2str(numBytes) ' bytes, ' num2str(numIter) ' Iterations. ' num2str(etime(tend,tstart)) ' seconds to sim.'];
title(ttl)
xlabel('E_{b}/N_{0} [dB]');
ylabel('Bit Error Rate [%]')
grid on
hold on
bert = berawgn(EbNodB,'psk',2,'nondiff');
plot(EbNodB,bert,'-*','LineWidth',2)
%plot(ebno0(1:9),bertCE(1:9),'-s','LineWidth',2)

legend('6 Mbps','9 Mbps','12 Mbps','18 Mbps','24 Mbps','36 Mbps','48 Mbps',...
    '54 Mbps','BPSK in awgn','Location','southwest')
%'BPSK with ecc in awgn','Location','southwest')