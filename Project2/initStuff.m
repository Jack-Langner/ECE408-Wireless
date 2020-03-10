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

numIter = 1e1;
tstart = clock;
wb = waitbar(0,'Beginning simulation');
for mm = 1:length(M)
for ii = 1:numIter
    for jj = 1:lenSNR
        newMsg1 = ['Diversity order: ' num2str(M(mm))];
        newMsg2 = ['Simulation iteration # ' num2str(ii)];
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

T = 2;
M = 2.^(0:1);
SNRdb2 = (0:40).';
lenSNR = length(SNRdb2);
%numIter = 1e2;
Aber = zeros(length(SNRdb2),length(M));
tstart = clock;
wb = waitbar(0,'Beginning simulation');
%set(wb,'Position',[10 10 10 10])

for mm = 1:length(M)
for ii = 1:numIter
    for jj = 1:lenSNR
        %SNR of ' num2str(SNRdb(jj)) 'dB,
        newMsg1 = ['Receive Diversity order: ' num2str(M(mm))];
        newMsg2 = ['Simulation iteration # ' num2str(ii)];
        %newMsg3 = 'how big is the box';
        newMsg = {newMsg1, newMsg2};
        waitbar( ((mm-1)*(lenSNR*numIter)+(ii-1)*lenSNR+jj)/(length(M)*numIter*lenSNR),wb,newMsg)
        if isequal(jj,1) && isequal(ii,4)
            tstart2 = clock;
        end

tmp = jackAlam(n,T,M(mm),SNRdb2(jj));
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

figure
semilogy(SNRdb,ber(:,1),'-o','LineWidth',2);
BERawgn = berawgn(SNRdb,'psk',2,'nondiff');
%semilogy(SNRdb,BERawgn,'LineWidth',3)
grid on
hold on
semilogy(SNRdb,ber(:,2),'-v','LineWidth',2);
% grid on
% hold on
semilogy(SNRdb,ber(:,3),'-s','LineWidth',2);

% last input for berfading is diversity order
%BERrf = NaN(length(SNRdb),1);

BERrf = berfading(SNRdb,'psk',2,1);
%semilogy(SNRdb,abs(BERrf),'LineWidth',3);
semilogy(SNRdb2,Aber(:,1),'-d','LineWidth',2)
semilogy(SNRdb2,Aber(:,2),'-^','LineWidth',2)
%semilogy(SNRdb,Aber(:,3),'-*','LineWidth',2)

%legend('Simulation','berfading','FontSize',24)
xlabel('SNR [dB]');ylabel('P_{b}, bit error rate (BER)');
%xlabel('SNR [dB]');ylabel('Bit Error Rate (BER) [%]');
legend('no diversity (1 Tx, 1 Rx)','MRRC (1 Tx, 2 Rx)',...
    'MRRC (1 Tx, 4 Rx)','new scheme (2 Tx, 1 Rx)',...
    'new scheme 2 Tx, 2 Rx)','FontSize',20)
% title('New Scheme and MRRC at equal power')
% legend('MRRC (1 Tx, 2 Rx)','MRRC (1 Tx, 4 Rx)','MRRC (1 Tx, 8 Rx)',...
%     'new scheme (2 Tx, 1 Rx)','new scheme (2 Tx, 2 Rx)',...
%     'new scheme (2 Tx, 4 Rx)','FontSize',20)
ylim([1e-6 1])
