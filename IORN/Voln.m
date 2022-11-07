
A=5;
a=1.5;
b=5;
x=-5:0.1:5;

Psi=A*exp(-x.^2/a^2+1i*b*x);
f=x;

Psi_rmlb=real(Psi);
Psi_rhnd=A*exp(-x.^2/a^2).*cos(x.*b);

Psi_2mlb=abs(Psi).^2;
Psi_2hnd=A^2*exp(-2*x.^2/a^2);

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
plot(x, Psi_rmlb, '-b', x, Psi_rhnd, '--or', x, Psi_2mlb, '-k', x, Psi_2hnd, '--gx');
xlabel('x values','Interpreter','latex');
ylabel('R($\psi$)','Interpreter','latex');
title('$\displaystyle R(\psi) \& {(\psi)}^2$','interpreter','latex');
legend('Psi_rmlb','Psi_rhnd','Psi_2mlb', 'Psi_2hnd')
grid on;

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
subplot(2, 2, 1), plot(x, Psi_rmlb, '-b');
xlabel('x values','Interpreter','latex');
ylabel('R($\psi$)','Interpreter','latex');
title('Действительная часть функции(mlb)');
grid on;
subplot(2, 2, 2), plot(x, Psi_rhnd,'--or');
xlabel('x values','Interpreter','latex');
ylabel('R($\psi$)','Interpreter','latex');
title('Действительная часть функции(hand)');
grid on;
subplot(2, 2, 3), plot(x, Psi_2mlb, '-k');
xlabel('x values','Interpreter','latex');
ylabel('R($\psi$)','Interpreter','latex');
title('Квадрат модуля волновой функции (mlb)');
grid on;
subplot(2, 2, 4), plot(x, Psi_2hnd, '--gx');
xlabel('x values','Interpreter','latex');
ylabel('R($\psi$)','Interpreter','latex');
title('Квадрат модуля волновой функции (hand)');
grid on;