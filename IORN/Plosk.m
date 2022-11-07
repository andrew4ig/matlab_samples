n=100;
lmbd=1;
c=3;
w=2*pi*c/lmbd;
t=-1:1/n:1;
x=-1:1/n:1;
k=2*pi/lmbd;
pfi=pi/4;
E0=1;

[X,T]=meshgrid(x,t);

E=E0*cos(T.*w-k.*X-pfi);
surf(X,T,E);
shading interp;
colorbar;
xlabel('Координата пространства x, м');
ylabel('Координата времени t, c');
zlabel('Плоская монохроматическая волна E=E(x,t)');

clear T E
T=0;
E=E0*cos(T*w-k*x-pfi);
figure
plot(x,E);
xlabel('Координата пространства x, м');
ylabel('E=E(x,0)');
title('Плоская монохроматическая волна E=E(x,0)');