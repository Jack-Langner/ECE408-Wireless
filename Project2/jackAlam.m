

function ber = jackAlam(n,T,M,SNRdb)

% For right now this function assumes BPSK modulated signals. n is the
% number of symbols to be transmitted, T the number of transmit antennas, M
% the number of the reeive antennas, SNRdB is self explanatory I hope.

% T = 3;
% M = 2;
% n = 5;
% SNRdb = 9;

b = randi(2,n,T)-1; %generate the data for number of antennas
c = -2*b+1; % BPSK modulate binary data
c = 2*c;
s = NaN(n*T,T);

% generate correct sequence of symbols given by Alamouti
s(1:T:end-1,:) = c;
s(2:T:end,:) = [-conj(c(:,2)) conj(c(:,1))];
%
% generate chan that is slowly changing, i.e. constant for two symbols
chan = genRayleighFadingV2(n,1,T*M,'false');
scc = NaN(n*T,T*M); %slowly changing channel
scc(1:T:end-1,:) = chan;
scc(2:T:end,:) = chan;

SCC = reshape(scc,n*T,T,M);
rx = NaN(n*T,M);
for mm = 1:M
    %tmp(:,mm) = sum(s.*SCC(:,:,mm),2);
    %s.*SCC is valid and works as expected
    rx(:,mm) = awgn(sum(s.*SCC(:,:,mm),2),SNRdb,'measured');
end
% %rx = awgn(tmp,SNRdb,'measured');%received data

st = NaN(n,T);

st(:,1) = sum(conj(chan(:,1:T:end-1)).*rx(1:T:end-1,:),2)...
    + sum(chan(:,2:T:end).*conj(rx(2:T:end,:)),2);

st(:,2) = sum(conj(chan(:,2:T:end)).*rx(1:T:end-1,:),2)...
    - sum(chan(:,1:T:end-1).*conj(rx(2:T:end,:)),2);

d = NaN(n,T,2);
d(:,:,1) = (st-1).*( conj(st)-conj(1) );
d(:,:,2) = (st-(-1)).*( conj(st)-conj((-1)) );

[~,ind] = min(d,[],3);
alam = ind-1;
ber = sum(abs(b-alam),'all')/(n*T); %'all' option beginning in 2018b

end