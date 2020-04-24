%%
t = [0:1e-3:pi/2 (pi/2)*ones(1,2e3)  pi/2:1e-3:pi];
w = sin(t).^2;
plot(1:length(t),w)
% 802.11 window function for OFDM
%% borrowing from myself
numSym = 280;
LFSR = ones(7,1);
pilotPolar = NaN(2^length(LFSR)-1,1);
for nn = 1:2^length(LFSR)-1
    nxt = mod(LFSR(end)+LFSR(4),2);
    pilotPolar(nn) = nxt;
    LFSR = [nxt; LFSR(1:end-1)];
end
clear LFSR
pilotPolar = -2*pilotPolar+1;
pilotPolarInd = mod(0:numSym-1,127)+1; % used for indexing purposes later

M = NaN(48,1);%NaN(rateStruct.NMSPOS,1);
LSI = (0:48-1).';%(0:rateStruct.NMSPOS-1).'; %logical subcarrier index
M(LSI>=0 & LSI<=4) = LSI(LSI>=0 & LSI<=4)-26;
M(LSI>=5 & LSI<=17) = LSI(LSI>=5 & LSI<=17)-25;
M(LSI>=18 & LSI<=23) = LSI(LSI>=18 & LSI<=23)-24;
M(LSI>=24 & LSI<=29) = LSI(LSI>=24 & LSI<=29)-23;
M(LSI>=30 & LSI<=42) = LSI(LSI>=30 & LSI<=42)-22;
M(LSI>=43 & LSI<=47) = LSI(LSI>=43 & LSI<=47)-21;
%
Pk = zeros(53,1);
pilotInds = [6 20 34 48];
pilotVals = [1 1 1 -1];
Pk(pilotInds) = pilotVals;
%% system parameters for 802.11n
dF = 20e6/64;
TFFT = 1/dF;
TG = TFFT/4;
dt = TFFT/64;
t = 0:dt:TFFT-dt;
k = (-26:26).';
k2 = (-32:31).';
r = NaN(64,1)+1j*NaN(64,1);

%% ifft without explicit for loop
d = (1:48).';
% for ii = 1:length(t)
% t1 = d.*exp(1j*2*pi*M*dF*(t(ii)-TG)); %first sum
% t2 = sum(t1);
% t3 = Pk.*exp(1j*2*pi*k*dF*(t(ii)-TG));
% t4 = pilotPolar(pilotPolarInd(1))*sum(t3);
% r(ii) = t2+t4;
% end
% new way of doing DFT without (explicit) for loops
T1 = sum(d.*exp(1j*2*pi*dF*M*(t-TG)));
T2 = pilotPolar(pilotPolarInd(1))*sum(Pk.*exp(1j*2*pi*dF*k*(t-TG)));
R = (T1+T2).'/8;

% subplot(1,2,1)
% plot(real(R))
% title('real(R)')
% subplot(1,2,2)
% plot(imag(R))
% title('imag(R)')

r = sum(R.*exp(-1j*2*pi*dF*k2*(t))).'/8;

%%

x = [0 0 0 0 0 0 1 -1 -1 1 -1 1 1 -1 -1 1 1 -1 1 -1 -1 -1 -1 -1 -1 1 -1 1 -1 1 -1 -1 ...
    0 1 -1 -1 -1 -1 -1 1 1 1 -1 -1 1 -1 -1 1 -1 -1 1 -1 -1 -1 1 -1 1 -1 -1 0 0 0 0 0].';
%% mandelbrot set attempt
n = 3;
d = 2;
dn = 1e-3;
z = (-n:dn:n)+1j*(n:-dn:-n).';
z1= z;
z2 = zeros([size(z) d]);
for ii = 1:d
z2(:,:,ii) = z1.^2+1;
%z22 = abs(z2)<=2;
z1 = z2(:,:,ii);
end
z22 = abs(sum(z2,3))<=2;
imagesc(z22)