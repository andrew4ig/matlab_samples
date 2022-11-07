datestr(now)
clear
clc
hbar=1973;
e=3.795;
m0=0.51e6;

Np=1000;
r=linspace(1e-10,10,Np);
dx=r(2)-r(1);

%H=U+Ep
U=zeros(1,Np);
for i=1:Np
    U(i)=-e^2/r(i);
end

koef=-hbar^2/(2*m0*12*dx^2);

E=eye(Np)*(-30);
for i=1:Np-1
    E(i,i+1)=E(i,i+1)+16;
    E(i+1,i)=E(i+1,i)+16;
end
for i=1:Np-2
    E(i,i+2)=E(i,i+2)-1;
    E(i+2,i)=E(i+2,i)-1;
end

E=E*koef;
%Hmiltonian
H=E+diag(U);

%finding eigenvalues and eigenvectors
[V,Ev]=eig(H);
E=diag(Ev);

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
plot(r, V(:,2), 'ko--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$?$', 'Interpreter', 'latex');

