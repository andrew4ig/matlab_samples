clear
% close all
clc
datestr(now)
load ('constants.mat','hbar','m0', 'J2eV');

%forming a task
e=1.6e-19;
L=101e-10;
Np=50;
dx=L/Np;
x=linspace(-L/2,L/2, Np);
Psi=@(n) sqrt(2/L)*sin((pi*n.*x)/L);
koef=-hbar^2/(2*m0*(dx^2));

%H=U+Ep
U=linspace(0,e,Np);
% U(1:Np/3)=0.01/J2eV;
% U(Np/3:2*Np/3)=0.1/J2eV;
% U=abs(x/(L/2)/J2eV/2);
% U(1)=e; U(Np)=0;

E=diag(ones(1,Np)*(-2))+diag(ones(1,Np-1),-1)+diag(ones(1,Np-1),1);
% E=eye(Np)*(-30);
% for i=1:Np-1
%     E(i,i+1)=E(i,i+1)+16;
%     E(i+1,i)=E(i+1,i)+16;
% end
% 
% for i=1:Np-2
%     E(i,i+2)=E(i,i+2)-1;
%     E(i+2,i)=E(i+2,i)-1;
% end

% E(Np,1)=16;
% E(1,Np)=16;
E(Np,1)=1;
E(1,Np)=1;

%Hamiltonian
H=E*koef+diag(U);
% E(1,Np)=1;
% E(Np,1)=1;

%finding eigenvalues and eigenvectors
[P,Eii]=eig(H);
Ei=diag(Eii);

N=1;

%Now we should notice that squared wave function has square that equals to dx
P=P*sqrt(1/dx);
% fprintf('Main state energy equals to %.2gJ=%.2g\n', [Ei(1), Ei(1)*J2eV]);
% fprintf('N=25 state energy equals to %.2gJ=%.2g\n', [Ei(N), Ei(N)*J2eV]);

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
subplot(2,2,1), 

plot(x, P(:,3)/max(abs(P(:,3))));
hold on;
plot(x, U/e, 'r--','LineWidth',2);
% legend('Numerical','Analitical');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$\Psi_1$', 'Interpreter', 'latex');
% xlim([0 L]);
grid on;

subplot(2,2,3), 
plot(x, P(:,N)/max(abs(P(:,N))));
hold on;
plot(x, U/e, '--c', 'LineWidth', 0.5)
%plot(x, Psi(N), 'r--','LineWidth',2);
% legend('Numerical','Analitical');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$\Psi_{25}$', 'Interpreter', 'latex');
% xlim([0 L]);
grid on;

subplot(2,2,2), 
plot(x, P(:,3).^2/max(abs(P(:,3).^2)));
hold on;
plot(x, U/e, '--c', 'LineWidth', 0.5)
%plot(x, Psi(1).^2, 'r--','LineWidth',2);
% legend('Numerical','Analitical');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$\Psi_1^2$', 'Interpreter', 'latex');
% xlim([0 L]);
grid on;

subplot(2,2,4), 
plot(x, P(:,N).^2/max(abs(P(:,N).^2)));
hold on;
plot(x, U/e, '--c', 'LineWidth', 0.5)
%plot(x, Psi(N).^2, 'r--','LineWidth',2);
% legend('Numerical','Analitical');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$\Psi_{25}^2$', 'Interpreter', 'latex');
% xlim([0 L]);
grid on;
