a1=[2 0 0];
a2=[1,2,0];
a3=[0,0,1];

la1=sqrt(sum(a1.*a1));
la2=sqrt(sum(a2.*a2));
la3=sqrt(sum(a3.*a3));
fprintf('la1=%g la2=%g la3=%g', [la1 la2 la3]);

V=sum(a1.*cross(a2,a3));
fprintf('V=%g',V);

b1=2*pi*cross(a2,a3)/V;
b2=2*pi*cross(a3,a1)/V;
b3=2*pi*cross(a1,a2)/V;
fprintf('b1=[%g %g %g]\n b2=[%g %g %g]\n b3=[%g %g %g]\n', [b1 b2 b3]);

lb1=sqrt(sum(b1.*b1));
lb2=sqrt(sum(b2.*b2));
lb3=sqrt(sum(b3.*b3));
fprintf('lb1=%g lb2=%g lb3=%g', [lb1 lb2 lb3]);

Vb=sum(b1.*cross(b2,b3));
fprintf('Vb=%g',Vb);

a1b1=sum(a1.*b1);
a1b2=sum(a1.*b2);
fprintf('a1b1=%g\na1b2=%g',[a1b1 a1b2]);

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
%hold on
plot3([0, a1(1), a2(1), a3(1)], ...
      [0, a1(2), a2(2), a3(2)], ...
      [0, a1(3), a2(3), a3(3)], ...
      'k.', 'LineWidth' , 7);
  
plot3([0, a1(1)],[0, a1(2)],[0, a1(3)], ...
      'k', 'LineWidth' , 2);
plot3([0, a2(1)],[0, a2(2)],[0, a2(3)], ...
      'k', 'LineWidth' , 2);
plot3([0, a3(1)],[0, a3(2)],[0, a3(3)], ...
      'k', 'LineWidth' , 2);

title('$hist(A)$', 'Interpreter', 'Latex');