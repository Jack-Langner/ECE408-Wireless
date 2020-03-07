%% ECE 408 - Wireless Communications
% Project 2 - Alamouti Transmit Diversity
% Jack Langner - MATLAB 2019b
% Due March 11, 2020

function r = genRayleighFading(N,fD,numChan)
r = NaN(N,numChan);
for qq = 1:numChan

% generate the samples of a Rayleigh fading channel with a given maximum
% Doppler frequency, fD. N is the number of points to be generated, r is a
% complex column vector that represents the channel at different times
%N = 1e1; %total number of noise samples

% N = 2^8;
% fD = 1;

ssbfdn = sqrt(1)*(randn(floor(N/2),2)+1j*randn(floor(N/2),2)); % single sideband frequency domain noise
%ssbfdn = (randn(floor(N/2),2));%+1j*randn(floor(N/2),2)); % single sideband frequency domain noise
if isequal(mod(N,2),0)
    fdn = [conj(flipud(ssbfdn));ssbfdn];
elseif isequal(mod(N,2),1)
    fdn = [conj(flipud(ssbfdn));zeros(1,2);ssbfdn];
end
%
% N = 2^8;
% fD = 1;
df = (2*fD)/(N-1);
T = 1/df;
f = (-fD:df:fD).';

SE = 1.5./(pi*fD*sqrt(1-(f/fD).^2)); %clark spectrum
dSE = (1.5*f)./(pi*(fD^3)*(1-(f/fD).^2).^(3/2));% derivative of clark spectrum

y0 = dSE(end-1)*(f(end)-f(end-1))+SE(end-1); %linear approx to SE(fD);
ind = [1 length(f)];
SE(ind) = [y0 y0];
sqSE = sqrt(SE);
% plot(f,SE)
% hold on
% plot(f,sqSE)

% first row of r will be real valued, because it is average.
q = fdn.*sqSE;
Q = ifft(q,N);
tmp = sqrt(sum(Q.^2,2));
rbar = mean(abs(tmp));
r(:,qq) = tmp./rbar;
end
% figure
% histogram(abs(r),50,'Normalization','pdf')
% x = 0:0.01:5;
% s2 = (mean(abs(r))/sqrt(pi/2))^2;
% y = (x/s2).*exp(-(x.^2)/(2*s2));
% hold on
% plot(x,y)
%plot(10*log10(abs(r)))

end