c=3e8; h=6.63e-34; k=1.38e-23; T=1500;
    lmbd=linspace(0.2e-6,6e-6);
    R=@(l)2*pi*c^2*h./l.^5./(exp(h*c./(l*k*T))-1);
    figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
    plot(lmbd,R(lmbd));
    ylabel('$R$','Interpreter', 'latex');
    xlabel('$\lambda, meters$','Interpreter', 'latex');
    title('$14. R={2pic^2h\over \lambda^5}{1\over e^{(h*c/(\lambda*k*T))}-1}$','Interpreter', 'latex');
    grid on;
    hold on;
    
    n=10^7;
    [l1, R1]=fminbnd(@(l)-R(l*n),0.2e-6*n,6e-6*n);
    plot(l1/n*10,R(l1/n*10),'ok', 'MarkerSize', 7,'MarkerFaceColor', 'b');
    text(l1/n*10,R(l1/n*10),['   Maxvalue',num2str(l1/n), '   ', num2str(-R1)]);
    
    [l1, R1]=fminbnd(@(l)-R(l),0.2e-6,6e-6);
    plot(l1,-R1,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'b');
    text(l1,-R1,['   Maxvalue',num2str(l1/n), '   ', num2str(-R1)]);
    
    RR=diff(R(lmbd))/(lmbd(2)-lmbd(1));
    [~ , ln]=min(RR);
    plot(lmbd(ln),R(lmbd(ln)),'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
    text(lmbd(ln),R(lmbd(ln)),'   Maxvalue');
    
    syms f(l)
    f(l)=R;
    Df=diff(f,l);
    Df=Df(lmbd);
    [~, l2]=min(Df);