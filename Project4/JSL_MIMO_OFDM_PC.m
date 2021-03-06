



function [ber1,ber2] = JSL_MIMO_OFDM_PC(rate,numBytes,snrdb,fD,CM)
% 2x2 MIMO with Precoding
% rate = 4;
% numBytes = 120;
% snrdb = 20;
% fD = 1;
nTX = 2;
nRX = 2;

[msgTX1, msgBin1] = JSL_OFDMModulate(rate,numBytes);
[msgTX2, msgBin2] = JSL_OFDMModulate(rate,numBytes);
nSamp = length(msgTX1);
x = [msgTX1 msgTX2].';

h = genRayleighFadingV3(nSamp,fD,nTX*nRX);%number of samples, fD, channels.

H = reshape(h.',nRX,nTX,[]);

U = NaN(nRX,nRX,nSamp)+1j*NaN(nRX,nRX,nSamp);
S = NaN(nRX,nTX,nSamp)+1j*NaN(nRX,nTX,nSamp);
V = NaN(nTX,nTX,nSamp)+1j*NaN(nTX,nTX,nSamp);
for ii = 1:nSamp
    [U(:,:,ii),S(:,:,ii),V(:,:,ii)] = svd(H(:,:,ii));%svd(eye(2));%
end

xt = NaN(nTX,nSamp);
for ii = 1:nSamp
    xt(:,ii) = V(:,:,ii)*x(:,ii);
end

%g = NaN(2,1,nSamp);
if isequal(CM,0)
    for ii = 1:nSamp
        g(:,:,ii) = H(:,:,ii)*xt(:,ii);
        %g(:,:,ii) = H*x(:,ii);
    end
elseif isequal(CM,1)
    g(:,:,1) = H(:,:,1)*xt(:,1);
    for ii = 2:nSamp
        g(:,:,ii) = H(:,:,ii)*xt(:,ii)+0.1*g(:,:,ii-1);
    end
elseif isequal(CM,2)
    g(:,:,1) = H(:,:,1)*xt(:,1);
    g(:,:,2) = H(:,:,2)*xt(:,2)+0.1*g(:,:,1);
    for ii = 3:nSamp
        g(:,:,ii) = H(:,:,ii)*xt(:,ii)+0.1*g(:,:,ii-1)+0.2*g(:,:,ii-2);
    end
end
sp = mean(abs(g).^2,3);
s2 = sp.*10.^(-snrdb/10);
n = sqrt(s2/2).*(randn(nRX,1,nSamp)+1j*randn(nRX,1,nSamp));

y = g+n;
yt = NaN(nRX,1,nSamp);
for ii = 1:nSamp
    yt(:,:,ii) = U(:,:,ii)'*y(:,:,ii);
end

y2 = reshape(yt(:),nRX,[]).';

msgDec1 = JSL_OFDMDemodulate(rate,numBytes,y2(:,1));
msgDec2 = JSL_OFDMDemodulate(rate,numBytes,y2(:,2));

ber1 = sum(abs(msgBin1-msgDec1))/numel(msgBin1);
ber2 = sum(abs(msgBin2-msgDec2))/numel(msgBin2);
%BER = [ber1 ber2];

end