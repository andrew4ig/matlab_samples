i=sqrt(-1);
tetha=linspace(0, pi);
phi=linspace(0,2*pi);
[TETHA, PHI]=meshgrid(tetha, phi);

R=sqrt(15/(32*pi)).*sin(TETHA)^3.*exp(i*3.*PHI);

[X, Y, Z]=sph2cart(PHI, pi/2-TETHA, real(R));
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
subplot(1, 3, 1), surf(X, Y, Z);
axis('equal');
colormap('gray');
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
zlabel('z','Interpreter','latex');
title('1', 'Interpreter','latex');
shading interp;

[X, Y, Z]=sph2cart(PHI, pi/2-TETHA, imag(R));
subplot(1, 3, 2), surf(X, Y, Z);
axis('equal');
colormap('gray');
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
zlabel('z','Interpreter','latex');
title('2', 'Interpreter','latex');
shading interp;

[X, Y, Z]=sph2cart(PHI, pi/2-TETHA, abs(R));
subplot(1, 3, 3), surf(X, Y, Z);
axis('equal');
colormap('gray');
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
zlabel('z','Interpreter','latex');
title('3', 'Interpreter','latex');
shading interp;