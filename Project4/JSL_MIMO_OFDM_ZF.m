

function BER = JSL_MIMO_OFDM_ZF(rate,numBytes,snrdb,fD)
% 2x2 MIMO with Zero Forcing (psuedo inverse)
% rate = 6;
% numBytes = 120;
% snrdb = 40;
% fD = 1;
nTX = 2;
nRX = 2;

[msgTX1, msgBin1] = JSL_OFDMModulate(rate,numBytes);
[msgTX2, msgBin2] = JSL_OFDMModulate(rate,numBytes);
nSamp = length(msgTX1);
x = [msgTX1 msgTX2].';

h = genRayleighFadingV3(nSamp,fD,nTX*nRX);%number of samples, fD, channels.
% w = [randperm(nSamp,nSamp);randperm(nSamp,nSamp);...
%      randperm(nSamp,nSamp);randperm(nSamp,nSamp)].';
% h = h0(w(:,1),:) 
H = reshape(h.',nRX,nTX,[]);
Wzf = NaN(nRX,nTX,nSamp); %only thing that doesnt scale with arbitrary nTX & nRX

for ii =1:nSamp
    %H(:,:,ii) = eye(2);
    Wzf(:,:,ii) = pinv(H(:,:,ii)); 
end
%
%H = eye(nRX,nTX)+sqrt(0.0001/2)*(randn(2)+1j*randn(2));
g = NaN(nRX,1,nSamp);
for ii = 1:nSamp
    g(:,:,ii) = H(:,:,ii)*x(:,ii);
    %g(:,:,ii) = H*x(:,ii);
end
sp = mean(abs(g).^2,3);
s2 = sp.*10.^(-snrdb/10);
n = sqrt(s2/2).*(randn(nRX,1,nSamp)+1j*randn(nRX,1,nSamp));

y = g+n;
yt = NaN(nRX,1,nSamp);
for ii = 1:nSamp
    yt(:,:,ii) = Wzf(:,:,ii)*y(:,:,ii);
end

y2 = reshape(yt(:),nRX,[]).';
msgDec1 = JSL_OFDMDemodulate(rate,numBytes,y2(:,1));
msgDec2 = JSL_OFDMDemodulate(rate,numBytes,y2(:,2));

ber1 = sum(abs(msgBin1-msgDec1))/numel(msgBin1);
ber2 = sum(abs(msgBin2-msgDec2))/numel(msgBin2);
BER = [ber1 ber2];
end