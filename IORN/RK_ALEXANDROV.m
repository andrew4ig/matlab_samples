clc
clear
close all
format shortG
%1
fprintf('1 task\ndate_n_time')
datestr (now)

%2

fprintf('2 task\n')
A=randi([2,6],3, 5);
disp('A=');
disp(A);
figure('Units', 'normalized','OuterPosition', [0 0 1/2 1])
subplot (1,2,1) 
histogram(reshape(A, 1, []))
title('$hist(A)$', 'Interpreter', 'Latex');

B=1+2*randn(5,3);
disp('B=');
disp(B);
subplot (1,2,2) 
histogram(reshape(B,1,[]))
title('$hist(B)$', 'Interpreter', 'Latex');

%3
fprintf('3 task\n')
[m, n]=size(A);
m=ceil(m/2);
n=ceil(n/2);

fprintf('\nA(%i,%i)=%g->',[m,n,A(m,n)])
A(m,n)=A(m,n)*(-0.5)-1;
fprintf('A(%i,%i)=%g',[m,n,A(m,n)])

[m, n]=size(B);
m=ceil(m/2);
n=ceil(n/2);

fprintf('\nB(%i,%i)=%g->',[m,n,B(m,n)])
B(m,n)=B(m,n)*(sind(45))+1;
fprintf('B(%i,%i)=%g\n',[m,n,B(m,n)])

%4
fprintf('\n4 task C=A.*B\n')
C=A.*B';
disp(C);

%5
fprintf('5 task D matrix \n')
[m, n]=size(C);
M=min(m, n);
N=max(m, n);
D=zeros(M);

for i=1:M
    D(i,i)=1;
end

for i=1:M-1
    D(i,i+1)=-0.5;
end
for i=1:M-1
    D(i+1,i)=sind(45);
end
disp('D=');
disp(D);

T=A*B+D;
disp('T=');
disp(T);

clear i
%6
fprintf('6 task - Basis\n')
%syms i j k;
a1=[2,0,0];
a2=[1,2,0];
a3=[0,0,1];

la1=sqrt(sum(a1.*a1));
la2=sqrt(sum(a2.*a2));
la3=sqrt(sum(a3.*a3));
fprintf('la1=%g la2=%g la3=%g\n', [la1 la2 la3]);
V=sum(a1.*cross(a2,a3));
fprintf('\nV=%g\n',V);
% format longEng
% disp(la1);
% disp(la2);
% disp(V);
% format short
% disp(la1);
% disp(la2);
% disp(V);

b1=2*pi*cross(a2,a3)/V;
b2=2*pi*cross(a3,a1)/V;
b3=2*pi*cross(a1,a2)/V;
fprintf('\nb1=[%g %g %g]\nb2=[%g %g %g]\nb3=[%g %g %g]\n', [b1 b2 b3]);

lb1=sqrt(sum(b1.*b1));
lb2=sqrt(sum(b2.*b2));
lb3=sqrt(sum(b3.*b3));
fprintf('\nlb1=%g lb2=%g lb3=%g', [lb1 lb2 lb3]);


Vb=sum(b1.*cross(b2,b3));
fprintf('\nVb=%g\n',Vb);
% format longEng
% disp(lb1);
% disp(lb2);
% disp(Vb);

a1b1=sum(a1.*b1);
a1b2=sum(a1.*b2);
fprintf('\na1b1=%g\na1b2=%g',[a1b1 a1b2]);
% format longEng
% disp(a1b1);
% disp(a1b2);

figure('Units', 'normalized', 'OuterPosition', [0 0 1/2 1]);
hold on
plot3([0, a1(1), a2(1), a3(1)], ...
      [0, a1(2), a2(2), a3(2)], ...
      [0, a1(3), a2(3), a3(3)], ...
      'k.', 'MarkerSize' , 7);
plot3([0, a1(1)],[0, a1(2)],[0, a1(3)], ...
      'k', 'LineWidth' , 2);
plot3([0, a2(1)],[0, a2(2)],[0, a2(3)], ...
      'k', 'LineWidth' , 2);
plot3([0, a3(1)],[0, a3(2)],[0, a3(3)], ...
      'k', 'LineWidth' , 2);
  
plot3([0, b1(1), b2(1), b3(1)], ...
      [0, b1(2), b2(2), b3(2)], ...
      [0, b1(3), b2(3), b3(3)], ...
      'ks', 'MarkerSize' , 7, ...
      'MarkerFaceColor', 'r');
plot3([0, b1(1)],[0, b1(2)],[0, b1(3)], ...
      'r--', 'LineWidth' , 2);
plot3([0, b2(1)],[0, b2(2)],[0, b2(3)], ...
      'r--', 'LineWidth' , 2);
plot3([0, b3(1)],[0, b3(2)],[0, b3(3)], ...
      'r--', 'LineWidth' , 2);
 
text(a1(1)+0.1,a1(2)+0.1,a1(3),'$\vec{a_1}$', 'Interpreter', 'latex');
text(a2(1)+0.1,a2(2)+0.1,a2(3),'$\vec{a_2}$', 'Interpreter', 'latex');
text(a3(1)+0.1,a3(2)+0.1,a3(3),'$\vec{a_3}$', 'Interpreter', 'latex');
text(b1(1)+0.1,b1(2)+0.1,b1(3),'$\vec{b_1}$', 'Interpreter', 'latex');
text(b2(1)+0.1,b2(2)+0.1,b2(3),'$\vec{b_2}$', 'Interpreter', 'latex');
text(b3(1)+0.1,b3(2)+0.1,b3(3),'$\vec{b_3}$', 'Interpreter', 'latex');
view(3);
title('$Structure$', 'Interpreter', 'Latex');
xlabel('The coordinate of Ox, i', 'Interpreter', 'latex');
ylabel('The coordinate of Oy, j', 'Interpreter', 'latex');
zlabel('The coordinate of Oz, k', 'Interpreter', 'latex');
grid on
axis equal




%7
fprintf('\n7 task - -Fermi-Dirak\n')
V=0.2;
temp=27;
t_to_T=@(t) t+273;
T=t_to_T(temp);
Ef=0;
e=1.6e-19;

kb_eV=8.617e-5;
ev_J=1.6e-19;
kb_J=kb_eV*ev_J;
E=linspace(-0.25, 0.25);

FD=@(Ef) 1./(1+exp((E-Ef)/(kb_eV*T)));
F=FD(E);
F1=FD(V/2*e/ev_J);
F2=FD(-V/2*e/ev_J);

figure('Units', 'normalized','OuterPosition', [0 0 1/2 1])
plot(F1, E, F2,E);
ylabel('$E(eV)$','Interpreter', 'Latex');
xlabel('$Fermi function f(E), T=27^{o}C$','Interpreter', 'Latex');
ylim([-0.25 0.25]);
hold on;
plot(abs(F1-F2), E, 'g--', 'LineWidth', 2);
legend('$V=0,1V$','$V=-0,1V$', '$|F_1-F_2|$', 'Interpreter', 'Latex');
grid on;

dx=E(2)-E(1);
F12=abs(F1-F2);
Fprods=F12.*dx;
Integ=sum(Fprods);
text(0.42, -0.22, ['S = ', num2str(Integ)], 'Interpreter', 'Latex');

