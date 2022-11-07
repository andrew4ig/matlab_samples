format short

o=3e-5;
n=250;
x=-o:o/n:o;
t=-o:o/n:o;
[X,T]=meshgrid(x,t);
i=sqrt(-1);

E=5*1.6e-19;
m=9.1e-31;
h=6.626e-34;
c=3e8;
l=c*h/E;
dl=l/100;
l=l-dl:dl/n:l+dl;
k=2*pi./l;
dk=abs(k(2*n+1)-k(1));
A=1;
w=2*pi*c./l;
v=sqrt(2*E/m);
vgroup=v;
vphaze=c^2/v;


T=0;
z=transpose(k);
F=A*cos(w.*T-transpose(k).*x);
F=sum(F)/dk;


figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
plot(x,F);
xlabel('x','Interpreter','latex');
ylabel('E(x,0)','Interpreter','latex');
xlim([-o, o])
text(x(100),550, ['\rightarrow v =', num2str(v)])
text(-2.5e-7,F(n+1)+10, ['\rightarrow v_g =', num2str(vgroup)])
text(x(n-85),F(n-85), ['\rightarrow v_p =', num2str(vphaze)])
text(x(n+85),F(n+85), ['\rightarrow v_g =', num2str(vphaze)])





