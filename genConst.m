
bpsk = [1 -1];

qpsk = ([-1 1] + 1j*[1;-1])/sqrt(2);
qpsk = qpsk(:);

qam16 = ((-3:2:3)+1j*(3:-2:-3).')/sqrt(10);
qam16 = qam16(:);

qam64 = ((-7:2:7)+1j*(7:-2:-7).')/sqrt(42);
qam64 = qam64(:);

qam256 = ((-15:2:15)+1j*(15:-2:-15).')/sqrt(170);
qam256 = qam256(:);
%%
figure
scatter(bpsk,[0 0],100,'x','LineWidth',2)
hold on
scatter(real(qpsk),imag(qpsk),100,'s','LineWidth',2)
scatter(real(qam16),imag(qam16),100,'h','LineWidth',2)
scatter(real(qam64),imag(qam64),100,'*')
%scatter(real(qam256),imag(qam256),'b')
legend('BPSK','QPSK','16QAM','64QAM')
grid on
xlabel('real axis');ylabel('imaginary axis')
title('802.11-2012 Modulation Constellations')