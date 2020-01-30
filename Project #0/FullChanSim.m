%% ECE 408 - Wireless Communications
% Project 0 - Channel Simulation
% Jack Langner - MATLAB 2018b
% Assigned January 22, 2020

%%
numIter = 10;
numSym = 1e3;
EbNodB = 0:0.5:12;
lenEbNo = length(EbNodB);

tPoly = poly2trellis([5 4],[23 35 0; 0 5 13]);
codeRate = 2/3;
traceBack = 12;
decDelay = 2*traceBack;

M = 2;
k = log2(M);
SNRdB = EbNodB+10*log10(k*codeRate);

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
        
        tBits = randi(2,numSym*k,1)-1;
        cBits = convenc(tBits,tPoly);
        cBitsMat = reshape(cBits,length(cBits)/k,k);

        tMsg = bi2de(cBitsMat);
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
           sqrt(10^(-SNRdB(jj)/10)/2)*(randn(numSym/codeRate,k)+1j*randn(numSym/codeRate,k));
        else
        rSig = awgn(txChan,SNRdB(jj),'measured');
        end
        
        if isequal(chan,1)
            rMsg = qamdemod(rSig,M);
        else
            dfeeq = dfe(8,4,lms(0.005),qammod(0:M-1,M)); %works best
            %dfeeq = dfe(8,4,varlms(0.005,1e-4,0.001,0.01),qammod(0:M-1,M));
            %dfeeq = dfe(8,4,normlms(0.005,0.01),qammod(0:M-1,M));
            %dfeeq = dfe(8,4,signlms(0.005),qammod(0:M-1,M));
            %dfeeq = dfe(10,5,rls(1),qammod(0:M-1,M));
            [symEst,yd] = equalize(dfeeq,rSig,tSig(1:TL));
            rMsg = qamdemod(yd,M);
        end
        rBitsMat = de2bi(rMsg);
        rcBits = reshape(rBitsMat,[],1);
        dBits = vitdec(rcBits,tPoly,traceBack,'cont','hard');
        [~, berVec(ii,jj)] = biterr(tBits(TL+1:end-decDelay),dBits(TL+decDelay+1:end));
   end
end
tend = clock;
close(wb)
fprintf('Time for simulation was %.2f seconds.\n',etime(tend,tstart))
ber = mean(berVec,1);
%
figure
semilogy(EbNodB,ber)
if isequal(M,2)
    bert = berawgn(EbNodB,'psk',2,'nondiff');
else
    bert = berawgn(EbNodB,'qam',M);
end
hold on
semilogy(EbNodB,bert,'r')
semilogy(qam4ecc,qam2ecc,'g')
if isequal(chan,1)
    ttl = ['Bit Error Rate for ' num2str(M) 'QAM, no ISI, ' num2str(numIter) ' Iterations. ' num2str(etime(tend,tstart)) ' seconds to sim.'];
else
    ttl = ['Bit Error Rate for ' num2str(M) 'QAM with ISI ' num2str(numIter) ' Iterations. ' num2str(etime(tend,tstart)) ' seconds to sim.'];
end
title(ttl)
legend('Simulated','Theoretical','Theorectical ECC','Location','southwest')
