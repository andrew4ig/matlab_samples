% This is a 1D FDTD simulation of psi. 

clear
clc

NN = 400;             % Number of points in the problem space.
hbar = 1.054e-34;     % Plank's constant
m0 = 9.1e-31;         % Free space mass of an electron
meff = 1.0;           % Effective mass: Si is 1.08, Ge is 0.067, GaAs is 0.55
melec = meff*m0;      % Mass of an electron
eV2J = 1.6e-19;       % Energy conversion factors
J2eV = 1/eV2J;

dx = .1e-9;           % The cell size
dt = 2e-17;           % Time steps 
ra = (0.5*hbar/melec)*(dt/dx^2);  % ra must be < .15
DX = dx*1e9;          % Cell size in nm.
XX = (DX:DX:DX*NN);   % Length in nm for plotting

% --- Specify the potential ---------------------
 
V=zeros(1,NN);
V(NN/2:NN)=0.1*eV2J;

% ---------------------------------------------------
% Initialize a sine wave in a gaussian envelope

lambda = 50;         % Pulse wavelength
sigma = 50;          % Pulse width
nc = 150;            % Starting position
prl = zeros(1,NN);   % The real part of the state variable
pim = zeros(1,NN);   % The imaginary part of the state variable
ptot = 0.;
for n=2:NN-1
    prl(n) = exp(-1.*((n-nc)/sigma)^2)*cos(2*pi*(n-nc)/lambda) ;
    pim(n) = exp(-1.*((n-nc)/sigma)^2)*sin(2*pi*(n-nc)/lambda) ;
    ptot = ptot +  prl(n)^2 + pim(n)^2;
end
pnorm = sqrt(ptot);    % Normalization constant

% Normalize and check
ptot = 0.;
for n=1:NN
    prl(n) = prl(n)/pnorm;
    pim(n) = pim(n)/pnorm;
    ptot = ptot +  prl(n)^2 + pim(n)^2;
end
ptot                   % This should have the value 1

T = 0;
n_step = 1;
while n_step > 0
    n_step = input('How many time steps  -->');

    % -----------This is the core FDTD program -------------
    for m=1:n_step
        T = T + 1;

        for n=2:NN-1
            prl(n) = prl(n) - ra*(pim(n-1) -2*pim(n) + pim(n+1)) ...
               + (dt/hbar)*V(n)*pim(n);
        end 

        for n=2:NN-1
            pim(n) = pim(n) + ra*(prl(n-1) -2*prl(n) + prl(n+1)) ...
               - (dt/hbar)*V(n)*prl(n);
        end

    end
    square=sum(prl.^2+pim.^2)
    prl=prl/sqrt(square);
    pim=pim/sqrt(square);
    square=sum(prl.^2+pim.^2)
    
    % ------------------------------------------------------

    % Calculate the expected values

    PE = 0.;
    for n=1:NN
        psi(n) = prl(n) + i*pim(n);    % Write as a complex function
        PE = PE + psi(n)*psi(n)'*V(n);
    end
    psi*psi'                % This checks normalization
    PE = PE*J2eV;           % Potential energy

    ke = 0. + j*0.;
    for n=2:NN-1
        lap_p = psi(n+1) - 2*psi(n) + psi(n-1);
        ke = ke + lap_p*psi(n)';
    end
    KE = -J2eV*((hbar/dx)^2/(2*melec))*real(ke); % Kinetic energy

    subplot(2,1,1)
    plot(XX,prl,'k') 
    hold on
    plot(XX,pim,'-.k')
    plot(XX,J2eV*V,'--k')
    hold off
    axis( [ 1 DX*NN -.2 .3 ])
    TT = text(5,.15,sprintf('%7.0f fs',T*dt*1e15));
    set(TT,'fontsize',12)
    TT = text(5,-.15,sprintf('KE = %5.3f eV',KE));
    set(TT,'fontsize',12)
    TT = text(25,-.15,sprintf('PE = %5.3f eV',PE));
    set(TT,'fontsize',12)
    TT = text(25,.13,sprintf('E_t_o_t = %5.3f eV',KE+PE));
    set(TT,'fontsize',12)
    xlabel('nm')
    set(gca,'fontsize',12)
    T

end
