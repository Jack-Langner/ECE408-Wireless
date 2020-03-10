%% ECE 408 - Wireless Communications
% Project 2 - Alamouti Transmit Diversity
% Jack Langner - MATLAB 2019b
% Due March 11, 2020

function r = genRayleighFadingV2(N,fD,numChan,fig)

% generate the samples of a Rayleigh fading channel with a given maximum
% Doppler frequency, fD. N is the number of points to be generated, r is a
% complex column vector that represents the channel at different times
%N = 1e1; %total number of noise samples
% improved speed by getting rid of for loop

ssbfdn = randn(floor(N/2),2*numChan)+1j*randn(floor(N/2),2*numChan); % single sideband frequency domain noise
% checking length of data
if isequal(mod(N,2),0)
    fdn = [conj(flipud(ssbfdn));ssbfdn];
elseif isequal(mod(N,2),1)
    fdn = [conj(flipud(ssbfdn));zeros(1,2);ssbfdn];
end
%
df = (2*fD)/(N-1);
f = (-fD:df:fD).'; % frquency vector

SE = 1.5./(pi*fD*sqrt(1-(f/fD).^2)); %Doppler spectrum
dSE = (1.5*f)./(pi*(fD^3)*(1-(f/fD).^2).^(3/2));
% derivative of Doppler spectrum

y0 = dSE(end-1)*(f(end)-f(end-1))+SE(end-1); %linear approx to SE(fD);
ind = [1 length(f)];
SE(ind) = [y0 y0];
sqSE = sqrt(SE);

% first row of r will be real valued, because it is average.
q = fdn.*sqSE;
Q = (ifft(q,N)).^2;
Q = reshape(Q,N,2,numChan);
Q = sqrt(reshape(sum(Q,2),N,numChan)); % the rayleigh noise


rbar = mean(abs(Q));
r = Q./rbar; % unit mean

if isequal(fig,'true')
figure % ploting a histogram vs pdf
histogram(abs(r(:,1)),50,'Normalization','pdf','LineWidth',2)
x = 0:0.01:3.5;
s2 = (mean(abs(r(:,1)))/sqrt(pi/2))^2;
y = (x/s2).*exp(-(x.^2)/(2*s2));
hold on
plot(x,y,'LineWidth',4)
legend('generated data','pdf','FontSize',24)
title('Rayleigh PDF')

figure % comparing the power spectrums
rPSD = abs(fft(r));
plot(f,SE,'LineWidth',2)
hold on
plot(f,mean(rPSD,2)/(max(mean(rPSD,2))/y0))
xlim(1.25*[-fD fD]);ylim([0 20])
xlabel('frequency [Hz]');ylabel('Magnitude');
title('Power Spectrum')
legend('Ideal','Generated','FontSize',24)
end

end