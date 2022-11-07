c=3e8; %Скорость света%
lam=transpose(linspace(525e-9,550e-9,500));%Длина рассматриваемой волны%
v=c./lam;%Вывод частоты%
k=2*pi./lam;%Волновое число%
A=1; %Амплитуда волны%
x=linspace(-300e-7,300e-7,500);
omega=1;
t=0;
S=A*exp(-1i*(omega*t-k.*x));
s1=real(S);
itog=sum(s1);
plot(x,itog);