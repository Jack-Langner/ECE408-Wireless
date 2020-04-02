%% ECE 408 - Wireless Communications
% Project 3 - Spread Spectrum, Secret Messages, and Cat Videos
% Jack Langner - MATLAB 2019b
% Due March whenever, 2020
% super rough draft
%%
rxdata = load('Rcvd_Langner.mat');
rxdata = rxdata.Rcvd; % received data vector

fml = 255; % frame length
usf = 4; % up sample factor / university of southern Florida

% correct implementation of m sequence with checks
n = 8;
initStates = de2bi((0:2^n-1).',n);
BPMS = NaN(2^n,2^n - 1); %BPSK modulated m sequence 
Mseq = BPMS;

for jj = 1:2^n
Tseq = NaN(1,2^n-1); % place holder
LFSR = initStates(jj,:);
for ii = 1:2^n-1
    nxt = mod(sum(LFSR([1 2 7 8])),2);
    Tseq(ii) = LFSR(1);
    LFSR = [LFSR(2:end) nxt];
end
BPMS(jj,:) = 1-2*Tseq;
Mseq(jj,:) = Tseq;
end
% [sum(Mseq,2) sum(~Mseq,2)] % check to see if valid M seq
% Rxcheck = BPMS*BPMS.'; % check autocorrelation of M seq
% imagesc(Rxcheck) % visualize the autocorrelation of M seq

clear LFSR nxt Tseq ii jj
%
usMS = zeros(2^n,fml*usf); %up sampled M seq
usMS(:,1:usf:end) = BPMS;
%
[B_RCOS, A_RCOS] = rcosine(1,4,'fir/sqrt',0.75);
%b = rcosdesign(0.75,6,4,'sqrt'); gives an equivalent filter

mfr = filter(B_RCOS,A_RCOS,rxdata);
%%
y = NaN(256,length(mfr));
for jj = 2:256
y(jj,:) = abs(filter(fliplr(usMS(jj,:)),1,mfr));
end

%plot(y)
%%

for jj = 129:129
    %jj
P = "";
    inds = find(y(jj,:)>240);
data = mfr(inds(1):inds(end)-1);

frames = length(data)/(fml*usf);
dsd = data.*repmat(usMS(jj,:),1,frames); %descambled data
ds2d = dsd(1:4:end); % down sampled, descrambled

ff = reshape(ds2d.',fml,[]);
h = hadamard(8);
for qq = 1:frames
t1 = ff(1:192,qq);
t2 = reshape(t1,8,[]);
pc = mean(angle(mean(t2)));
t3 = t2*exp(-1j*pc);

d1 = t3.'*h(:,6);
%
d2 = zeros(size(d1));
d2(real(d1)<0)=1;
d3 = (reshape(d2,8,[])).';
%d4 = bi2de(d3,'left-msb');
d4 = bi2de(d3,'right-msb');
p = convertCharsToStrings(char(d4));
P = strcat(P,p);
end
P
end