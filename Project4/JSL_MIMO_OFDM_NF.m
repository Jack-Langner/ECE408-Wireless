

function BER = JSL_MIMO_OFDM_NF(rate,numBytes,snrdb,fD)
% 2x2 MIMO no channel filter
% rate = 6;
% numBytes = 150;
% snrdb = 40;
% fD = 1;
nTX = 2;
nRX = 2;

[msgTX1, msgBin1] = JSL_OFDMModulate(rate,numBytes);
[msgTX2, msgBin2] = JSL_OFDMModulate(rate,numBytes);
nSamp = length(msgTX1);
x = [msgTX1 msgTX2].';

%
h = genRayleighFadingV3(nSamp,fD,nTX*nRX);%number of samples, fD, channels.
% h = h(w); 
H = reshape(h.',2,2,[]);
%H = eye(nRX,nTX)+sqrt(0.0001/2)*(randn(2)+1j*randn(2));
g = NaN(2,1,nSamp);
for ii = 1:nSamp
    g(:,:,ii) = H(:,:,ii)*x(:,ii);
    %g(:,:,ii) = H*x(:,ii);
end
sp = mean(abs(g).^2,3);
s2 = sp.*10.^(-snrdb/10);
n = sqrt(s2/2).*(randn(2,1,nSamp)+1j*randn(2,1,nSamp));

y = g+n;
y2 = reshape(y(:),2,[]).';
msgDec1 = JSL_OFDMDemodulate(rate,numBytes,y2(:,1));
msgDec2 = JSL_OFDMDemodulate(rate,numBytes,y2(:,2));

ber1 = sum(abs(msgBin1-msgDec1))/numel(msgBin1);
ber2 = sum(abs(msgBin2-msgDec2))/numel(msgBin2);
BER = [ber1 ber2];
end