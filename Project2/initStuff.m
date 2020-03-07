%% ECE 408 - Wireless Communications
% Project 2 - Alamouti Transmit Diversity
% Jack Langner - MATLAB 2019b
% Due March 11, 2020
%%

M = 2.^(0:2);
n = 1e4;
SNRdb = (0:50).';
s2 = 10.^(-SNRdb/10);
lenSNR = length(SNRdb);
ber = zeros(length(SNRdb),length(M));

numIter = 10;
tstart = clock;
wb = waitbar(0,'Beginning simulation');
%set(wb,'Position',[10 10 10 10])
for mm = 1:length(M)
for ii = 1:numIter
    for jj = 1:lenSNR
        %SNR of ' num2str(SNRdb(jj)) 'dB,
        newMsg1 = ['Diversity order: ' num2str(M(mm))];
        newMsg2 = ['Simulation iteration # ' num2str(ii)];
        %newMsg3 = 'how big is the box';
        newMsg = {newMsg1, newMsg2};
        waitbar( ((mm-1)*(lenSNR*numIter)+(ii-1)*lenSNR+jj)/(length(M)*numIter*lenSNR),wb,newMsg)
        if isequal(jj,1) && isequal(ii,4)
            tstart2 = clock;
        end

tmp = jackMRRC(n,M(mm),SNRdb(jj));
ber(jj,mm) = ber(jj,mm)+tmp;
if isequal(jj,1) && isequal(ii,4)
    tend2 = clock;
    fprintf('Estimated time for sim is %f seconds.\n',etime(tend2,tstart2)*numIter*lenSNR)
end
    end
end
end
tend = clock;
close(wb)
fprintf('Time for RX div simulation was %.2f seconds.\n',etime(tend,tstart))
ber = ber/numIter;
%%
%************ alamouti test************************************************
n = 1e4;
T = 2;
M = 2.^(0:1);
SNRdb = (0:50).';
lenSNR = length(SNRdb);
numIter = 50;
Aber = zeros(length(SNRdb),length(M));
tstart = clock;
wb = waitbar(0,'Beginning simulation');
%set(wb,'Position',[10 10 10 10])

for mm = 1:length(M)
for ii = 1:numIter
    for jj = 1:lenSNR
        %SNR of ' num2str(SNRdb(jj)) 'dB,
        newMsg1 = ['Diversity order: ' num2str(M(mm))];
        newMsg2 = ['Simulation iteration # ' num2str(ii)];
        %newMsg3 = 'how big is the box';
        newMsg = {newMsg1, newMsg2};
        waitbar( ((mm-1)*(lenSNR*numIter)+(ii-1)*lenSNR+jj)/(length(M)*numIter*lenSNR),wb,newMsg)
        if isequal(jj,1) && isequal(ii,4)
            tstart2 = clock;
        end

tmp = jackAlam(n,T,M(mm),SNRdb(jj));
Aber(jj,mm) = Aber(jj,mm)+tmp;
if isequal(jj,1) && isequal(ii,4)
    tend2 = clock;
    fprintf('Estimated time for sim is %f seconds.\n',etime(tend2,tstart2)*numIter*lenSNR)
end
    end
end
end
tend = clock;
close(wb)
fprintf('Time for TX div simulation was %.2f seconds.\n',etime(tend,tstart))
Aber = Aber/numIter;
%%
%**************** plotting ***********************************************


%
figure
semilogy(SNRdb,ber);

bert = berawgn(SNRdb,'psk',2,'nondiff');
grid on
hold on
%semilogy(SNRdb,bert)

%s1 = 10.^(-SNRdb/20);
qBERT = (qfunc(sqrt(2*10.^(SNRdb/10))));%theoretical bit error rate for BPSK
%semilogy(SNRdb,qBERT)

% use bertool for theoretical curves
% last input for berfading is diversity order
BERrf = NaN(length(SNRdb),1);
BERrf(:,1) = berfading(SNRdb,'psk',2,1);
% BERrf(:,1) = berfading(SNRdb,'psk',2,1);
% BERrf(:,2) = berfading(SNRdb,'psk',2,2);
% BERrf(:,3) = berfading(SNRdb,'psk',2,4);
semilogy(SNRdb,abs(BERrf));
semilogy(SNRdb,Aber)
xlabel('SNR [dB]');ylabel('BER [%]');
legend('sim 1','sim 2','sim 4','Rayleigh','Alam1','Alam2')
ylim([1e-6 1])
%legend('sim','berawgn','qfunc','Rayleigh')