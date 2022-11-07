clear
clc
datestr(now)
format shortg
load ('constants.mat','hbar','m0', 'J2eV');

%forming a task
L=40e-9;
dx=0.1e-9;
Np=L/dx;
x=linspace(0,L, Np);

dt=0.02e-15;
Mp=3400;

U=zeros(1, 400);

koefx=hbar*dt/(2*m0*(dx^2));
koeft=dt/hbar;

m=[1 1700 3400];
ms=length(m);
m=sort(m);

% ms = input('How many moments u want to see?\nms=');
% clear m
% fprintf('What are they?\n');
% for ms=1:ms
%     m(ms)=input(['m(', num2str(ms), ')=']);
% end
% m=sort(m);

A=400;
b=0.79e-9;
a=3.5e-9;
x0=15e-9;

G=A*exp(-((x-x0)/(a)).^2+1i*(x-x0)/b);
Gre=real(G);
Gim=imag(G);
G=G/sqrt(sum(Gre.^2+Gim.^2));
Gre=real(G);
Gim=imag(G);

figure('Units','Normalized', 'OuterPosition',  [0 0 1 1]);
k=1;

for Time=1:m(ms)
    for n=2:Np-1
       Gre(n)=Gre(n)-koefx*(Gim(n-1)-2*Gim(n)+Gim(n+1))+koeft*U(n)*Gim(n);
    end 
    
    for n=2:Np-1
       Gim(n)=Gim(n)+koefx*(Gre(n-1)-2*Gre(n)+Gre(n+1))-koeft*U(n)*Gre(n);
    end
    
    G2=Gre.^2+Gim.^2;
    if m(k)== Time
        subplot (ms,1,k) 
        plot(x, Gre,'r', x, Gim, 'b--', x, G2, 'k', x, U, 'c');
        xlabel('Coordinate X, m','Interpreter','latex');
        ylabel('{Psi}','Interpreter','latex');
        legend('Real','Imaginary','Density of probability', 'U');
        text(x(100),max(Gre)/2, ['t=', num2str(Time*0.02), ' femptosec'],'Interpreter','latex')
        k=k+1;
    end
end

