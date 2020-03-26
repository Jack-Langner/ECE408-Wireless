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
%frxd = reshape(rxdata.',fml*usf,[]); %formatted received data

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
% [sum(Mseq) sum(~Mseq)] % check to see if valid M seq
% Rxcheck = BPMS*BPMS.'; % check autocorrelation of M seq
% imagesc(Rxcheck) % visualize the autocorrelation of M seq

clear initStates LFSR nxt Tseq ii jj

usMS = zeros(2^n,fml*usf); %up sampled M seq
usMS(:,1:usf:end) = BPMS;
%
[B_RCOS, A_RCOS] = rcosine(1,4,'fir/sqrt',0.75);

mfr = filter(B_RCOS, A_RCOS,rxdata); %matched filter response
fmfr = reshape(mfr.',fml*usf,[]); %formatted mfr
numFrame = length(mfr)/(fml*usf);

c1 = NaN(2^n,usf);

    for ii = 1:usf
        %c1(jj,ii) = abs(mfr(1:fml*4)*circshift((usMS(jj,:)).',ii-1));
        %c1(jj,ii) = abs(circshift(usMS(jj,:),ii-1)*fmfr(:,1));
        c1(:,ii) = abs(circshift(usMS,ii-1,2)*fmfr(:,1));
    end
clear ii
[max1 , ind1] = max(c1);

% plot(c1)
% legend('0','1','2','3')

[~,iind1] = max(max1);

%z = rxdata(4:4:fml*4);

% % shifting the data frame, probably not good
% c2 = NaN(2^8,10);
% for jj = 1:2^8
% for ii = 1:10
%     c2(jj,ii) = abs(mfr(ii:fml*4+ii-1)*(usMS(jj,:))');
% end
% end
% clear ii jj
% 
% [max2, ind2] = max(c2);
% metaD2 = [max2;ind2];
%%
%%
n = 8;
h = hadamard(n);
b = -(h-1)/2;
t1 = 0:n-1;
t2 = bi2de(de2bi(t1,'left-msb'))+1;
o = h(t2,:);
%%
% from Timmy H
%fml = 255;
SYM_RATE_TX = 1e6;%1e6/8;
SAMP_RATE_TX = 4e6;%500e3;
[B_RCOS, A_RCOS] = rcosine(SYM_RATE_TX,SAMP_RATE_TX,'fir/sqrt',0.75);

%mfr = conv(qqq, B_RCOS); % matched filter response
%rxds = (mfr(1:4:end-length(B_RCOS)+1)).'; %received down sampled signal

mfr = filter(B_RCOS, A_RCOS,rxdata);
rxds = (mfr(1:4:end)).'; % received, down sampled
%
M = (reshape(rxds,fml,length(rxds)/fml)).';
%initFrame = M(1,:);
%%
w = abs(initFrame*BPMS.');
[~, ind] = max(w);
%
%z = initFrame.*chck(ind,:)
z = BPMS(ind,:);
%y = initFrame.*chck(ind,:);

Y = M.*z;

ii = 1;
if2 = Y(ii,1:192); % intermediate frame 
pc = mean(unwrap(angle(if2))); %phase correction

if3 = if2*exp(-1j*pc); % trying to correct phase

q = ones(1,length(if3));
q(real(if3)<0)=-1; % thresholding
q = (reshape(q.',8,[])).';
%
TW = q*h; % test the walsh spectrum
imagesc(TW)
%x = [initFrame;z;y];
%%

tmp = abs(M*(BPMS.'));
[mx, inds] = max(tmp.');