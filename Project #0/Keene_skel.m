

numIter = 2e2;
numSym = 1e3;
SNRdb = 0:16;
lenSNR = length(SNRdb);

k = 2;
M =2^k;

chan = 1;
% chan = [1 0.2 0.4];
% chan = [0.227 0.460 0.688 0.460 0.227];

% ts = 1e-3;
% chan = rayleighchan(ts,1);
% chan.pathdelays = ts*(0:2);
% chan.AvgPathGaindB = 0:5:10;
% chan.StoreHistory = 1; %Uncomment if you want to be able to do plot(chan)

berVec = NaN(numIter,lenSNR);

wb = waitbar(0,'Performing simulation');
tstart = clock;
for ii = 1:numIter
    for jj = 1:lenSNR
        
        waitbar((lenSNR*(ii-1)+jj)/(numIter*lenSNR),wb)
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
        
        %noise = randn
        rSig = awgn(txChan,SNRdb(jj),'measured');
        rMsg = qamdemod(rSig,M);
        rBits = de2bi(rMsg);
        
        [~,berVec(ii,jj)] = biterr(tBits,rBits);
    end
end
tend = clock;
close(wb)
fprintf('Time for simulation was %f seconds.\n',etime(tend,tstart))

ber = mean(berVec,1);
semilogy(SNRdb,ber)

bert = berawgn(SNRdb,'psk',2,'nondiff');
hold on
semilogy(SNRdb+3.01,bert,'r')
legend('BER','BERT')