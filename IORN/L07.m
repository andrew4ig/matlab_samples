clear
clc
datestr(now)
load ('constants.mat','hbar','m0', 'J2eV', 'kT');

%forming a task
L=101e-10;
meff=0.067;
m0=m0*meff;
Np=100;
dx=L/Np;
x=linspace(0,L, Np);
Psi=@(n) sqrt(2/L)*sin((pi*n.*x)/L);
En = @(n) (pi*hbar)^2/(2*m0*(L)^2)*(n^2);
deltaEnergy = @(n) (pi*hbar)^2/(2*m0*(L)^2)*((n+1).^2-n.^2);
koef=-hbar^2/(2*m0*(dx^2));

%H=U+Ep
U=zeros(1,Np);

E=eye(Np)*(-2);
for i=1:Np-1
    E(i,i+1)=E(i,i+1)+1;
    E(i+1,i)=E(i+1,i)+1;
end
% 
% for i=1:Np-2
%     E(i,i+2)=E(i,i+2)-1;
%     E(i+2,i)=E(i+2,i)-1;
% end

%Hamiltonian
H=E*koef+diag(U);

%finding eigenvalues and eigenvectors
[P,Eii]=eig(H);
Ei=diag(Eii);

%Now we should notice that squared wave function has square that equals to dx
P=P*sqrt(1/dx);
fprintf('\nt_0 equals to %.2g\n', [koef*J2eV]);
fprintf('Main state energy numerical equals to %.2g\n', [Ei(1)*J2eV]);
fprintf('N=2 state energy numerical equals to %.2g\n', [Ei(2)*J2eV]);
fprintf('N=3 state energy numerical equals to %.2g\n', [Ei(3)*J2eV]);
fprintf('N=25 state energy numerical equals to %.2g\n', [Ei(25)*J2eV]);

fprintf('\nMain state energy analitical equals to %.2g\n', [En(1)*J2eV]);
fprintf('N=2 state energy analitical equals to %.2g\n', [En(2)*J2eV]);
fprintf('N=3 state energy analitical equals to %.2g\n', [En(3)*J2eV]);
fprintf('N=25 state energy analitical equals to %.2g\n', [En(25)*J2eV]);

fprintf('\nMain state energy delta equals to %.2g\n', [deltaEnergy(1)*J2eV*1000]);
fprintf('N=2 state energy delta equals to %.2g\n', [deltaEnergy(2)*J2eV*1000]);
fprintf('N=3 state energy delta equals to %.2g\n', [deltaEnergy(3)*J2eV*1000]);
fprintf('N=25 state energy delta equals to %.2g\n', [deltaEnergy(25)*J2eV*1000]);

fprintf('\n\n');

fprintf('kT equals to %.2g\n', [kT*1000]);
fprintf('\nt_0 equals to %.2g\n', [koef*J2eV]);
M=[ Ei(1)*J2eV Ei(2)*J2eV Ei(3)*J2eV Ei(25)*J2eV;
    En(1)*J2eV En(2)*J2eV En(3)*J2eV En(25)*J2eV;
    deltaEnergy(1)*J2eV*1000 deltaEnergy(2)*J2eV*1000 deltaEnergy(3)*J2eV*1000 deltaEnergy(25)*J2eV*1000]



fprintf('Maximum value of Psi_1_squared equals to %.2g\n', [max(P(:,1).^2)]);
