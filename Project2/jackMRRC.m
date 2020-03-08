

function ber = jackMRRC(n,M,SNRdb)
% For right now this function assumes BPSK modulated data transmitted from
% a single antenna. The data goes through a Rayleigh fading channel and
% white gaussian noise is added at the receiver. 
% n is the block size
% M is the number of receive antennas
% SNRdb is the signal to noise ratio in dB

% SNRdb = 20;
% M = 1;
% n = 1e3;

b = randi(2,n,1)-1;
s = -2*b+1;
chan = genRayleighFadingV2(n,1,M,'false');

rxData = NaN(n,M);
for qq = 1:M
    rxData(:,qq) = awgn(s.*chan(:,qq),SNRdb,'measured');
end
% SNRL = 10.^(-SNRdb/10);
% rx2 = s.*chan;
% mrx = mean(rx2);
% nois = sqrt(SNRL*mrx).*randn(n,M);
% rxData = rx2+nois;

st = sum(conj(chan).*rxData,2);

d = NaN(n,2);
d(:,1) = (st-1).*(conj(st)-conj(1));
d(:,2) = (st-(-1)).*(conj(st)-conj((-1)));

[~,ind] = min(d,[],2);
mrrc = ind-1;
ber = sum(abs(b-mrrc))/n;

end