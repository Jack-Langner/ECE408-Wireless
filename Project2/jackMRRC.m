

function ber = jackMRRC(n,M,SNRdb)
% For right now this function assumes BPSK modulated data transmitted from
% a single antenna. The data goes through a Rayleigh fading channel and
% white gaussian noise is added at the receiver. 
% n is the block size
% M is the number of receive antennas
% SNRdb is the signal to noise ratio in dB

b = randi(2,n,1)-1; %randomly generated data
s = -2*b+1; % BPSK modulate
chan = genRayleighFadingV2(n,1,M,'false'); % rayleigh channel

rxData = NaN(n,M);
for qq = 1:M
    rxData(:,qq) = awgn(s.*chan(:,qq),SNRdb,'measured');
end

st = sum(conj(chan).*rxData,2); %generate the s tilde's

d = NaN(n,2);
d(:,1) = (st-1).*(conj(st)-conj(1));
d(:,2) = (st-(-1)).*(conj(st)-conj((-1))); 
%calculating distance from s tilde to BPSK symbols

[~,ind] = min(d,[],2); %decision rule
mrrc = ind-1; 
ber = sum(abs(b-mrrc))/n; %

end