clear
clc
datestr(now)
load ('constants.mat','hbar','m0', 'J2eV');

%forming a task
L=101e-10;
Np=5000;
dx=L/Np;
x=linspace(0,L, Np);
Psi=@(n) sqrt(2/L)*sin((pi*n.*x)/L);
koef=-hbar^2/(2*m0*12*(dx^2));

%H=U+Ep
U=zeros(1,Np);

E=eye(Np)*(-30);
for i=1:Np-1
    E(i,i+1)=E(i,i+1)+16;
    E(i+1,i)=E(i+1,i)+16;
end

for i=1:Np-2
    E(i,i+2)=E(i,i+2)-1;
    E(i+2,i)=E(i+2,i)-1;
end

% E(Np,1)=16;
% E(1,Np)=16;
% E(Np-1,1)=-1;
% E(1,Np-1)=-1;

%Hamiltonian
H=E*koef+diag(U);

%finding eigenvalues and eigenvectors
[P,Eii]=eig(H);
Ei=diag(Eii);

N=25;

%Now we should notice that squared wave function has square that equals to dx
P=P*sqrt(1/dx);
fprintf('Main state energy equals to %.2gJ=%.2g\n', [Ei(1), Ei(1)*J2eV]);
fprintf('N=25 state energy equals to %.2gJ=%.2g\n', [Ei(N), Ei(N)*J2eV]);

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
subplot(2,2,1), 
plot(x, abs(P(:,2))/P(:,2)P(:,1));
hold on;
plot(x, Psi(1), 'r--','LineWidth',2);
legend('Numerical','Analitical');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$\Psi_1$', 'Interpreter', 'latex');
xlim([0 L]);
grid on;

subplot(2,2,3), 
plot(x, P(:,N));
hold on;
plot(x, Psi(N), 'r--','LineWidth',2);
legend('Numerical','Analitical');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$\Psi_{25}$', 'Interpreter', 'latex');
xlim([0 L]);
grid on;

subplot(2,2,2), 
plot(x, P(:,1).^2);
hold on;
plot(x, Psi(1).^2, 'r--','LineWidth',2);
legend('Numerical','Analitical');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$\Psi_1^2$', 'Interpreter', 'latex');
xlim([0 L]);
grid on;

subplot(2,2,4), 
plot(x, P(:,N).^2);
hold on;
plot(x, Psi(N).^2, 'r--','LineWidth',2);
legend('Numerical','Analitical');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$\Psi_{25}^2$', 'Interpreter', 'latex');
xlim([0 L]);
grid on;
