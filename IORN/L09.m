clear
clc
datestr(now)
format shortG
load ('constants.mat','hbar','m0', 'J2eV', 'eV2J');

%forming a task
L=40e-9;
dx=0.1e-9;
Np=L/dx;
x=linspace(0,L, Np);

dt=0.02e-15;
Mp=3400;

U=zeros(1, Np);
U(1:Np/2)=0.1*eV2J;

koefx=hbar*dt/(2*m0*(dx^2));
koeft=dt/hbar;

m=[1 1700 3400 3000];
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
b=3e-9/2/pi;
%a=120*dx/6;
a=3e-9;
x0=10e-9;

G=A*exp(-((x-x0)/(a)).^2+1i*(x-x0)/b);
Gre=real(G);
Gim=imag(G);
G=G/sqrt(sum(Gre.^2+Gim.^2)*dx);
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
    Gre=Gre/sqrt(sum(G2)*dx);
    Gim=Gim/sqrt(sum(G2)*dx);
    G2=Gre.^2+Gim.^2;
    
    EP = 0;
    GP = Gre + 1i.*Gim;       % Write as a complex function
    GP = GP/sqrt(GP*GP');
    for n=1:Np
        EP = EP + GP(n)*GP(n)'.*U(n);
    end
    EP = EP*J2eV;                % Potential energy
    
    EK = 0;
    for n=2:Np-1
        LAP = GP(n+1) - 2*GP(n) + GP(n-1);
        EK = EK + LAP*GP(n)';
    end
    EK = -J2eV*((hbar/dx)^2/(2*m0))*real(EK); % Kinetic energy

    if m(k)== Time
        subplot (ms,1,k) 
        plot(x, Gre,'r', x, Gim, 'b--',x, G2*5e-5, 'k', x, U*J2eV*1e5,'--g');
        xlabel('Coordinate X, m',   'Interpreter','latex');
        ylabel('{Psi}',             'Interpreter','latex');
        text(0.5e-9,1.5e4, sprintf('t=%5.2f fs',   Time*dt*1e15), 'Interpreter','latex')
        text(0.5e-9,0.5e4, sprintf('EK=%5.3f meV', EK*1000),      'Interpreter','latex')
        text(0.5e-9,-.5e4, sprintf('EP=%5.3f meV', EP*1000),      'Interpreter','latex')
        text(0.5e-9,-1.5e4,sprintf('E=%5.3f meV',  (EK+EP)*1000), 'Interpreter','latex')
        legend('Real','Imaginary','Density of probability', 'Force field','Interpreter','latex');
        ylim([-2e4,2e4]);
        xlim([0, 40e-9]);
        k=k+1;
    end
end

v=10e-9/68e-15;
Ek=(m0*v^2)/2;
Ek=Ek*J2eV*1000;

lmbd=5e-9;
p=2*pi*hbar/lmbd;
E = p^2/(2*m0)*J2eV*1000;