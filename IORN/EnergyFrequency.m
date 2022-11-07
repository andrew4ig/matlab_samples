function Z = EnergyFrequency(meff, L, n)
    load ('constants.mat', 'hbar', 'm0', 'J2eV', 'eV2J');
    En_mev=(pi*n*hbar).^2/(2*m0*meff*(L*1e-9)^2);
    Z=En_mev;
    wn=En_mev/hbar;
    disp(['For an electron meff=',num2str(meff),'in L=', num2str(L),'nm:']);
    fprintf('E%1i=%.4g meV; w%1i=%.3g rad/s\n', [n; En_mev*J2eV*1000; n; wn]);  
end

