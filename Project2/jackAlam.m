

function ber = jackAlam(n,T,M,SNRdb)
% For right now this function assumes BPSK modulated signals. n is the
% number of symbols to be transmitted, T the number of transmit antennas, M
% the number of the reeive antennas, SNRdB is self explanatory I hope.
% T = 2;
% M = 1;
% n = 5;
% SNRdb = 19;

b = randi(2,n,T)-1;
c = -2*b+1;
s = NaN(n*T,T);

s(1:T:end-1,:) = c;
s(2:T:end,:) = [-conj(c(:,2)) conj(c(:,1))];

chan = genRayleighFading(n,1,T*M);
c2 = NaN(n*T,T*M);
c2(1:T:end-1,:) = chan;
c2(2:T:end,:) = chan;

%rxData = NaN(n*T,M);

tmp = sum(s.*c2,2);

rxData = awgn(tmp,SNRdb,'measured');

st = NaN(n,T);
st(:,1) = conj(c2(1:T:end-1,1)).*rxData(1:T:end-1)+c2(1:T:end-1,2).*conj(rxData(2:T:end));

st(:,2) = conj(c2(1:T:end-1,2)).*rxData(1:T:end-1)-c2(1:T:end-1,1).*conj(rxData(2:T:end));

d = NaN(n,T,2);
d(:,:,1) = (st-1).*(conj(st)-conj(1));
d(:,:,2) = (st-(-1)).*(conj(st)-conj((-1)));

ind = NaN(n,T);
%alam= NaN(n,T);
for tt = 1:T
    [~,ind(:,tt)] = min(d(:,:,tt),[],2);
end
alam = ind-1;

ber = sum(abs(b-alam),'all');
%[alam b]
