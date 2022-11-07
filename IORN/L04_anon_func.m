load ('constants.mat', 'hbar', 'm0', 'J2eV', 'kb');
deltaEnergy = @(meff, L, n) (pi*hbar)^2/(2*m0*meff*(L*1e-9)^2)*((n+1).^2-n.^2);
Lmax=@(meff, T) sqrt((3*(pi*hbar)^2)/(2*kb*T*m0*meff));

MEFF=0.07;
L=20;
N=1:5;
T=300;

answ1 = deltaEnergy (MEFF, L, N);
disp(['For an electron meff=',num2str(MEFF),'in L=', num2str(L),'nm:'])
fprintf('deltaE%1i%1i = %.2g meV;\n', [N; N+1; answ1*J2eV*1000]);

answ2 = Lmax(MEFF, T);
fprintf(['\nFor practical system: \nmeff = ', num2str(MEFF),', T=', num2str(T),'; n=min(n)=1\n']);
fprintf('Lmax = %.2g nm;\n', [answ2*1e9]);


datestr(now)