%  syms x y
x=0:0.1:700;
y=(1.4/138.9).*((x/138.9).^(1.4-1)).*exp(-(x/138.9).^1.4);
plot(x,y);