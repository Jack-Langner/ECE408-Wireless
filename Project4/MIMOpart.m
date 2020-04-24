%% ECE 408 - Wireless Communications
% Project 4 - MIMO OFDM
% MIMO part
% Jack Langner - MATLAB 2019b
% Due April 29, 2020

%% System Parameters
dF = 20e6/64;
TFFT = 1/dF;
TG = TFFT/4;
dt = TFFT/64;
t = 0:dt:TFFT-dt;
k = (-26:26).';
k2 = (-32:31).';
r = NaN(64,1)+1j*NaN(64,1);

%% generate nakagami rand vars.
% in Statistics and Machine Learning Toolbox
pd = makedist('Nakagami','mu',5,'omega',2);
r = random(pd,1e4,1);
mean(r)
histogram(r,100,'Normalization','pdf')

%% QAM16
M = 16;
n = sqrt(M);
const = (-n:2:n)+1j*(n:-2:-n).'; 
acp = sum(abs(const).^2,'all')/numel(const);
constV = const(:)/sqrt(acp);%unit power
rd = constV(randi(16,48,1)); %randomly generated data
msp = mean(abs(rd).^2); %mean symbol power

%% 16PSK
MO = 16;
m = 0:MO-1;
const = exp(1j*2*pi*m/MO);
acp = sum(abs(const).^2)/numel(const);

plot(const,'*')
hold on
plot(cos(0:0.001:2*pi),sin(0:0.001:2*pi),'r')