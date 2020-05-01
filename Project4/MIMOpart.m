%% ECE 408 - Wireless Communications
% Project 4 - MIMO OFDM
% MIMO part
% Jack Langner - MATLAB 2019b
% Due April 29, 2020

%% System Parameters
dF = 20e6/64;
TFFT = 1/dF;
TG = TFFT/4;
dt = TFFT/64;
t = 0:dt:TFFT-dt;
k = (-26:26).';
k2 = (-32:31).';
r = NaN(64,1)+1j*NaN(64,1);

%% 2x2 mimo
% nTx = 2;
% nRx = 2;
% nSamp = 10;
% csm = eye(nTx,nRx)+sqrt(0.01/2)*(randn(nTx,nRx)+1j*randn(nTx,nRx));
% % channel state matrix
% [U,S,V] = svd(csm);
% TXD = randi(3,nSamp,nTx);
% 
% 
% RXD = TXD*csm;
%% 2x2 MIMO with 'normal' notation
snrdb = 20;
nTX = 2;
nRX = 2;
nSamp = 10;

d = randi(2,nTX,nSamp)-1; % data in bits
x = 1-2*d; %BPSK modulated
H = eye(nRX,nTX)+sqrt(0.01/2)*(randn(nRX,nTX)+1j*randn(nRX,nTX));
n = sqrt(0.01/2)*(randn(nRX,nSamp)+1j*randn(nRX,nSamp));
y = H*x+n;
%% 2x2 MIMO no channel filter
snrdb = 40;
nTX = 2;
nRX = 2;
nSamp = 1e3;
MO = 2;
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

h = genRayleighFadingV3(nSamp,1,4);%number of samples, fD, channels.
% w = [randperm(nSamp,nSamp);randperm(nSamp,nSamp);...
%      randperm(nSamp,nSamp);randperm(nSamp,nSamp)].';
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
check1 = var(n,0,3);

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
BER = sum(abs(binC-binRX),'all')/numel(binC)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2x2 MIMO with Precoding
snrdb = 20;
nTX = 2;
nRX = 2;
nSamp = 1e1;
MO = 256;
k = log(MO)/log(2);
p = sqrt(MO);
n = 2.^(0:k-1).';
%numSym = 5000;
%binData = randi(2,numSym*48,k)-1;
binData = randi(2,nSamp,k,nTX)-1;
decData = NaN(nSamp,nTX);
for ii = 1:nTX
decData(:,ii) = (binData(:,:,ii)*n)+1;
%decData(:,2) = (binData(:,:,2)*n)+1;
end
binC = dec2bin(decData(:)-1)-'0';
%decData = reshape(decData,nSamp,[]);
const = (-(p-1):2:(p-1))+1j*((p-1):-2:-(p-1)).'; 
%acp = sum(abs(const).^2,'all')/numel(const);
constV = const(:);%/sqrt(acp);%unit power
x = constV(decData).';

h = genRayleighFadingV3(nSamp,1,nTX*nRX);%number of samples, fD, channels.
% w = [randperm(nSamp,nSamp);randperm(nSamp,nSamp);...
%      randperm(nSamp,nSamp);randperm(nSamp,nSamp)].';
% h = h0(w(:,1),:) 
H = reshape(h.',nRX,nTX,[]);

U = NaN(nRX,nRX,nSamp)+1j*NaN(nRX,nRX,nSamp);
S = NaN(nRX,nTX,nSamp)+1j*NaN(nRX,nTX,nSamp);
V = NaN(nTX,nTX,nSamp)+1j*NaN(nTX,nTX,nSamp);
for ii = 1:nSamp
    [U(:,:,ii),S(:,:,ii),V(:,:,ii)] = svd(H(:,:,ii));
end

xt = NaN(nTX,nSamp);
for ii = 1:nSamp
    xt(:,ii) = V(:,:,ii)*x(:,ii);
end

%
%H = eye(nRX,nTX)+sqrt(0.0001/2)*(randn(2)+1j*randn(2));
g = NaN(nRX,1,nSamp);
for ii = 1:nSamp
    g(:,:,ii) = H(:,:,ii)*xt(:,ii);
    %g(:,:,ii) = H*x(:,ii);
end
sp = mean(abs(g).^2,3);
s2 = sp.*10.^(-snrdb/10);
n = sqrt(s2/2).*(randn(nRX,1,nSamp)+1j*randn(nRX,1,nSamp));
check1 = var(n,0,3);

y = g+n;
yt = NaN(nRX,1,nSamp);
for ii = 1:nSamp
    yt(:,:,ii) = U(:,:,ii)'*y(:,:,ii);
end

y2 = reshape(yt(:),nRX,[]);

dm = NaN(nRX,nSamp,length(constV));
for ii = 1:length(constV)
    dm(:,:,ii) = abs(y2-constV(ii));
end
[~,ind] = (min(dm,[],3));
i2 = ind.';

binRX = dec2bin(i2(:)-1)-'0';
%
BER = sum(abs(binC-binRX),'all')/numel(binC)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2x2 MIMO with Zero Forcing (psuedo inverse)
snrdb = 20;
nTX = 2;
nRX = 2;
nSamp = 1e3;
MO = 4;
k = log(MO)/log(2);
p = sqrt(MO);
n = 2.^(0:k-1).';
%numSym = 5000;
%binData = randi(2,numSym*48,k)-1;
binData = randi(2,nSamp,k,nTX)-1;
decData = NaN(nSamp,nTX);
for ii = 1:nTX
decData(:,ii) = (binData(:,:,ii)*n)+1;
%decData(:,2) = (binData(:,:,2)*n)+1;
end
binC = dec2bin(decData(:)-1)-'0';
%decData = reshape(decData,nSamp,[]);
const = (-(p-1):2:(p-1))+1j*((p-1):-2:-(p-1)).'; 
%acp = sum(abs(const).^2,'all')/numel(const);
constV = const(:);%/sqrt(acp);%unit power
x = constV(decData).';

h = genRayleighFadingV3(nSamp,1,nTX*nRX);%number of samples, fD, channels.
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
check1 = var(n,0,3);

y = g+n;
%
yt = NaN(nRX,1,nSamp);
for ii = 1:nSamp
    yt(:,:,ii) = Wzf(:,:,ii)*y(:,:,ii);
end

y2 = reshape(yt(:),nRX,[]);
%
dm = NaN(nRX,nSamp,length(constV));
for ii = 1:length(constV)
    dm(:,:,ii) = abs(y2-constV(ii));
end
[~,ind] = (min(dm,[],3));
i2 = ind.';

binRX = dec2bin(i2(:)-1)-'0';
%
BER = sum(abs(binC-binRX),'all')/numel(binC)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2x2 MIMO with Minimum Mean Square Error (MMSE) 
snrdb = 20;
nTX = 2;
nRX = 2;
nSamp = 1e3;
MO = 16;
k = log(MO)/log(2);
p = sqrt(MO);
n = 2.^(0:k-1).';
%numSym = 5000;
%binData = randi(2,numSym*48,k)-1;
binData = randi(2,nSamp,k,nTX)-1;
decData = NaN(nSamp,nTX);
for ii = 1:nTX
decData(:,ii) = (binData(:,:,ii)*n)+1;
%decData(:,2) = (binData(:,:,2)*n)+1;
end
binC = dec2bin(decData(:)-1)-'0';
%decData = reshape(decData,nSamp,[]);
const = (-(p-1):2:(p-1))+1j*((p-1):-2:-(p-1)).'; 
%acp = sum(abs(const).^2,'all')/numel(const);
constV = const(:);%/sqrt(acp);%unit power
x = constV(decData).';

h = genRayleighFadingV3(nSamp,1,nTX*nRX);%number of samples, fD, channels.
% w = [randperm(nSamp,nSamp);randperm(nSamp,nSamp);...
%      randperm(nSamp,nSamp);randperm(nSamp,nSamp)].';
% h = h0(w(:,1),:) 
H = reshape(h.',nRX,nTX,[]);
Wmmse = NaN(nRX,nTX,nSamp); %only thing that doesnt scale with arbitrary nTX & nRX

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
check1 = var(n,0,3);
for ii =1:nSamp
    %H(:,:,ii) = eye(2);
    Wmmse(:,:,ii) = (H(:,:,ii)'*H(:,:,ii)+diag(s2))\H(:,:,ii)';
end
y = g+n;
%
yt = NaN(nRX,1,nSamp);
for ii = 1:nSamp
    yt(:,:,ii) = Wmmse(:,:,ii)*y(:,:,ii);
end

y2 = reshape(yt(:),nRX,[]);
%
dm = NaN(nRX,nSamp,length(constV));
for ii = 1:length(constV)
    dm(:,:,ii) = abs(y2-constV(ii));
end
[~,ind] = (min(dm,[],3));
i2 = ind.';

binRX = dec2bin(i2(:)-1)-'0';
%
BER = sum(abs(binC-binRX),'all')/numel(binC)
%% channel filters ??
[U,S,V] = svd(H); %Precoding
xt = V*x;
yt = U'*(H*xt+n);

Wzf = pinv(H); %Zero Forcing Filter

Wmmse = ((H'*H+0.01*eye(nTX))^(-1))*H';
%%
rayChan = comm.RayleighChannel;
ricChan = comm.RicianChannel;
x = ones(1e4,1);
y1 = step(rayChan,x);
y2 = step(ricChan,x);

histogram(abs(y2),100,'Normalization','pdf')
%% generate nakagami rand vars.
% in Statistics and Machine Learning Toolbox
% pd = makedist('Nakagami','mu',5,'omega',2);
% r = random(pd,1e4,1);
% mean(r)
% histogram(r,100,'Normalization','pdf')
% 
% m = (gamma(5+1/2)/gamma(5))*sqrt(2/5);
% 
% R = fft(r,2^16);
% plot(abs(R))

%% QAM16
M = 16;
n = sqrt(M);
const = (-n:2:n)+1j*(n:-2:-n).'; 
acp = sum(abs(const).^2,'all')/numel(const);
constV = const(:)/sqrt(acp);%unit power
rd = constV(randi(16,48,1)); %randomly generated data
msp = mean(abs(rd).^2); %mean symbol power

%% 16PSK
MO = 16;
m = 0:MO-1;
const = exp(1j*2*pi*m/MO);
acp = sum(abs(const).^2)/numel(const);

plot(const,'*')
hold on
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'r')