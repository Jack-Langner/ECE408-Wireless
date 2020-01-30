%% ECE 408 - Wireless Communications
% Project 0 - Channel Simulation
% Jack Langner - MATLAB 2018b
% Assigned January 22, 2020

%%
numIter = 30;
numSym = 1e3;
EbNodB = 0:0.5:12;
lenEbNo = length(EbNodB);

M = 4;
k = log2(M);
SNRdB = EbNodB+10*log10(k);
%chan = 1;
chan = [1 0.2 0.4];
% chan = [0.227 0.460 0.688 0.460 0.227];

if isequal(chan,1)
    TL = 0;
else
    TL = 200;
end
berVec = NaN(numIter,lenEbNo);

wb = waitbar(0,'Performing simulation');
tstart = clock;
for ii = 1:numIter
    for jj = 1:lenEbNo
        
        waitbar((lenEbNo*(ii-1)+jj)/(numIter*lenEbNo),wb)
        %tBits = [ones(TL,k); randi(2,numSym-TL,k)-1];
        tBits = randi(2,numSym,k)-1;
        tMsg = bi2de(tBits);
        tSig = qammod(tMsg,M);
        
        if isequal(chan,1)
            txChan = tSig;
        elseif isa(chan,'channel.rayleigh')
            reset(chan);
            txchan = filter(chan,tSig);
        else
            txChan = filter(chan,1,tSig);
        end
        
        if isequal(M,2)
        rSig = txChan + ...
           sqrt(10^(-SNRdB(jj)/10)/2)*(randn(numSym,k)+1j*randn(numSym,k));
        else
        rSig = awgn(txChan,SNRdB(jj),'measured');
        end
        
        if isequal(chan,1)
            rMsg = qamdemod(rSig,M);
        else
            dfeeq = dfe(8,4,lms(0.005),qammod(0:M-1,M));
            [symEst,yd] = equalize(dfeeq,rSig,tSig(1:TL));
            rMsg = qamdemod(yd,M);
        end
        
        rBits = de2bi(rMsg);
        
        [~,berVec(ii,jj)] = biterr(tBits(TL+1:end,:),rBits(TL+1:end,:));
    end
end
tend = clock;
close(wb)
fprintf('Time for simulation was %.2f seconds.\n',etime(tend,tstart))
ber = mean(berVec,1);
%
figure
semilogy(EbNodB,ber,'-*')
if isequal(M,2)
    bert = berawgn(EbNodB,'psk',2,'nondiff');
else
    bert = berawgn(EbNodB,'qam',M);
end
hold on
semilogy(EbNodB,bert,'r')

if isequal(chan,1)
    ttl = ['Bit Error Rate for ' num2str(M) 'QAM, no ISI, ' num2str(numIter) ' Iterations. ' num2str(etime(tend,tstart)) ' seconds to sim.'];
else
    ttl = ['Bit Error Rate for ' num2str(M) 'QAM with ISI ' num2str(numIter) ' Iterations. ' num2str(etime(tend,tstart)) ' seconds to sim.'];
end
title(ttl)
legend('Simulated','Theoretical','Location','southwest')
