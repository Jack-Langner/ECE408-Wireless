
n = 1e4;
s2 = 3.32;
v = sqrt(s2)*randn(n,2);

r = sqrt(sum(v.^2,2));

histogram(r,50,'Normalization','pdf')

x = 0:0.01:10;

y = (x/s2).*exp(-(x.^2)/(2*s2));
hold on
plot(x,y)

rbar = mean(r);