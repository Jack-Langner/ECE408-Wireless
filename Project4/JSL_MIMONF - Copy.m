

function BER = JSL_MIMONF(nSamp,MO,snrdb,fD,CM)
% 2x2 MIMO no channel filter
% nSamp = 1e3;
% MO = 16;
% snrdb = 40;
% fD = 1;
% CM = 1;
nTX = 2;
nRX = 2;

if isequal(MO,2)
    constV = [1;-1];
    decData = randi(2,nSamp,nTX);
    binData = decData-1;
    binC = binData(:);
else
    k = log(MO)/log(2);
    p = sqrt(MO);
    n = 2.^(0:k-1).';
    const = (-(p-1):2:(p-1))+1j*((p-1):-2:-(p-1)).';
    constV = const(:);%/sqrt(acp);%unit power
    binData = randi(2,nSamp,k,nTX)-1;
    decData = NaN(nSamp,nTX);
    decData(:,1) = (binData(:,:,1)*n)+1;
    decData(:,2) = (binData(:,:,2)*n)+1;
    binC = dec2bin(decData(:)-1)-'0';
end

x = constV(decData).';
%
h = genRayleighFadingV3(nSamp,fD,nTX*nRX);%number of samples, fD, channels.
% w = [randperm(nSamp,nSamp);randperm(nSamp,nSamp);...
%      randperm(nSamp,nSamp);randperm(nSamp,nSamp)].';
% h = h(w);
H = reshape(h.',2,2,[]);
%H = eye(nRX,nTX)+sqrt(0.0001/2)*(randn(2)+1j*randn(2));
%CM = 2;
g = NaN(2,1,nSamp);
if isequal(CM,0)
    for ii = 1:nSamp
        g(:,:,ii) = H(:,:,ii)*x(:,ii);
        %g(:,:,ii) = H*x(:,ii);
    end
elseif isequal(CM,1)
    g(:,:,1) = H(:,:,1)*x(:,1);
    for ii = 2:nSamp
        g(:,:,ii) = H(:,:,ii)*x(:,ii)+0.1*g(:,:,ii-1);
    end
elseif isequal(CM,2)
    g(:,:,1) = H(:,:,1)*x(:,1);
    g(:,:,2) = H(:,:,2)*x(:,2)+0.1*g(:,:,1);
    for ii = 3:nSamp
        g(:,:,ii) = H(:,:,ii)*x(:,ii)+0.1*g(:,:,ii-1)+0.2*g(:,:,ii-2);
    end
end

sp = mean(abs(g).^2,3);
s2 = sp.*10.^(-snrdb/10);
n = sqrt(s2/2).*(randn(2,1,nSamp)+1j*randn(2,1,nSamp));
check1 = var(n,0,3);


%
y = g+n;
y2 = reshape(y(:),2,[]);
%
dm = NaN(nRX,nSamp,length(constV));
for ii = 1:length(constV)
    dm(:,:,ii) = abs(y2-constV(ii));
end
[~,ind] = (min(dm,[],3));
i2 = ind.';
binRX = dec2bin(i2(:)-1)-'0';
%
BER = sum(abs(binC-binRX),'all')/numel(binC);
end