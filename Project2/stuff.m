

N = 2^8;
fD = 1;
df = (2*fD)/(N-1);
f = (-fD:df:fD).';

SEz = (1.5)./(pi*fD*sqrt(1-(f/fD).^2));

dSEz= ((1.5)/(pi*fD))*(f./((fD^2)*((1-(f/fD).^2).^(3/2))));

y0 = dSEz(end-1)*df + SEz(end-1);
ti= [1 length(f)]; %temporary indices
SEz(ti) = [y0 y0];
sqSEz = sqrt(SEz);
%plot(f,[SEz; sqSEz])

cgn = randn(N/2,2)+1j*randn(N/2,2);

qqq = [conj(fliplr(cgn)); cgn];

sp = qqq.*sqSEz;
tdn = ifft(sp,N);

r = sqrt(sum(tdn.^2,2));

%%
n = 2;
T = 2;
b = -2*(randi(2,n,T)-1)+1 + 1j*(-2*(randi(2,n,T)-1)+1)

s = NaN(n*T,T);
s(1:T:end-1,:) = b;
s(2:T:end,:) = [-conj(b(:,2)) conj(b(:,1))];
%%
tmp = randi(5,3,1)+1j*randi(5,3,1)

TMP = toeplitz(tmp)

%%

s = sym('s',[3 1]);
h = sym('h',[3 1]);
n = sym('n',[3 1]);

r = toeplitz(-s)*h+n;

st = toeplitz(-h)*r

