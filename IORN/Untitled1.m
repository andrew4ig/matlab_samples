format short

o=3e-5;
n=100;
x=-o:o/n:o;
t=-o:o/n:o;
[X,T]=meshgrid(x,t);
i=sqrt(-1);

E=5*1.6e-19;
m=9.1e-31;
h=6.626e-34;
c=3e8;
l=c*h/E;
k=2*pi/l;
dk=k/100;
A=1;
v=sqrt(2*E/m);
vgroup=v;
vphaze=c^2/v;

k=linspace(k-dk,k+dk,2*n+1);
dk=k(2)-k(1);
F=E*exp(-i*(E*2*pi./h.*T-k.*X))*dk;
F=sum(F, 1);
F=real(F);



figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
plot(x,F);
xlabel('x','Interpreter','latex');
ylabel('E(x,0)','Interpreter','latex');
xlim([-o, o])
text(x(100),550, ['\rightarrow v =', num2str(v)])
text(-2.5e-7,F(n+1)+10, ['\rightarrow v_g =', num2str(vgroup)])
text(x(n-85),F(n-85), ['\rightarrow v_p =', num2str(vphaze)])
text(x(n+85),F(n+85), ['\rightarrow v_g =', num2str(vphaze)])