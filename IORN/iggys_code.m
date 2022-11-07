clear
close all
clc


k = 8.617e-5;   % постоянная Больцмана (эВ/К)
m0=9.10956e-31; % масса покоя электрона (кг)
e=1.60219e-19;  % заряд электрона (Кл)
hbar=1.054e-34; % постоянная Дирака (Дж с)
K2C=273.15;
meff=0.1*m0;

syms E U0 kl hb m a n positive;
lambda=0.5;
E=lambda*U0;
kl=sqrt((2*m*E)/hb^2);
%из условия сшивки получим%
l=simplify(kl*a==pi*n-2*asin((hb*kl)/sqrt(2*m*U0)));
z=solve(l,U0);
%6ерем значение от pi/2 до pi%
z=z*a^2;
z=subs(z,n,1);
z=vpa(subs(z,[hb,m],[hbar,meff]))*1e18/e;
fprintf('При lambda=%1.4g: a2U= %1.4g eV*nm^2',lambda,z);


l=simplify(kl*a==pi*n-2*asin((hb*kl)/sqrt(2*m*U0)));
N=solve(l,n);
