

x = (-5:0.01:5)';
p = -4:0.5:4;
y = sin(2*pi*(x-p))./(2*pi*(x-p));

plot(x,abs(y))
hold on
plot(x,abs(sum(y,2)),'--')

%%
n = 1e1;
x = randi(2,n,1)-1;
v = -2*x+1;

y = -(v-1)/2;
%%
%str = 'Hello Timmy H!';
str = 'Joy, bright spark of divinity, Daughter of Elysium, Fire-inspired we tread';
q = (double(str))';
a = dec2hex(q);
Q = de2bi(q,8);
W = bi2de(Q);
w = char(W);
ww = '';
for ii = 1:length(w)
    ww = [ww w(ii)];
end

%%

x = (0:63)';
v = ifft(x,64);
y = fft(v,64);
%%

i = 7:-2:-7;
q = (7:-2:-7)';
c = i+1j*q;
c = c(:);
mean(abs(c).^2)
%%
