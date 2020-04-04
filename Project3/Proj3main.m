%% ECE 408 - Wireless Communications
% Project 3 - Spread Spectrum, Secret Messages, and Cat Videos
% Jack Langner - MATLAB 2019b
% Due April 1, 2020, whoops
%%
rxdata = load('Rcvd_Langner.mat');
rxdata = rxdata.Rcvd; % received data vector

FL = 255; % frame length
usf = 4; % up sample factor / university of southern Florida
n = 8;
initStates = de2bi((0:2^n-1).',n);
BPMS = NaN(2^n,2^n - 1); %BPSK modulated m sequence

for jj = 1:2^n
    Tseq = NaN(1,2^n-1); % place holder
    LFSR = initStates(jj,:);
    for ii = 1:2^n-1
        nxt = mod(sum(LFSR([1 2 7 8])),2);
        Tseq(ii) = LFSR(1);
        LFSR = [LFSR(2:end) nxt];
    end
    BPMS(jj,:) = 1-2*Tseq; % BPSK modulate the M sequence
end
% [sum(Mseq,2) sum(~Mseq,2)] % check to see if valid M seq
Rxcheck = BPMS*BPMS.'; % check autocorrelation of M seq
imagesc(Rxcheck) % visualize the autocorrelation of M seq
title('PN Sequence Correlation')
clear initStates LFSR nxt Tseq ii jj
%
usMS = zeros(2^n,FL*usf); %up sampled M seq
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
%%
plot(y(129,:))
hold on
plot([1 length(y(129,:))],[240 240],'--k')
ylim([0 300])
title('PN and Data Correlation')
[~,ms] = max(y(:,1044)); % look in column 1044 based on advice from Timmy H
%plot(y)
%%
vc = [32 65:127]; %valid characters
P = "";
inds = find(y(ms,:)>240);
data = mfr(inds(1):inds(end)-1);

frames = length(data)/(FL*usf);
dsd = data.*repmat(usMS(ms,:),1,frames); %descambled data
ds2d = dsd(1:4:end); % down sampled, descrambled

ff = reshape(ds2d.',FL,[]);
fff = ff(1:192,:);

h = hadamard(8);
pc = NaN(frames,1);
for qq = 1:frames
    t1 = fff(:,qq);
    t2 = reshape(t1,8,[]);
    pc(qq) = mean(angle(mean(t2)));
    t3 = t2*exp(-1j*pc(qq));
    
    d1 = t3.'*h(:,6);
    d2 = zeros(size(d1));
    d2(real(d1)<0)=1;
    d3 = (reshape(d2,8,[])).';
    d4 = bi2de(d3,'right-msb');
    p = convertCharsToStrings(char(d4));
    P = strcat(P,p);
end
P
%%
pc2 = unwrap(pc);
chipRate = 1e6;
FR = chipRate/FL;
dpc2 = pc2(2:end)-pc2(1:end-1);
tm = (1/FR:1/FR:34/FR)*1e3;


plot(tm.',[pc2 pc].','-o','LineWidth',3)
title('Average Phase Frame to Frame')
legend('unwrapped','wrapped','Location','northwest','FontSize',20)
xlabel('time [ms]');ylabel('\theta [radians]')

figure
plot(tm(1:end-1),dpc2,'-o','LineWidth',3)
title('Derivative of Average Phase Frame to Frame')
xlabel('time [ms]');
ylabel('$$\frac{d\theta}{dt}$$ [radians/ms]','Interpreter','latex')

freq = (dpc2(1)*FR)/(2*pi)  