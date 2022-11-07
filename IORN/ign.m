p=100;
ax=1e-10;
ay=1e-10;
x=0:ax/p:ax;
y=0:ay/p:ay;
[X,Y]=meshgrid(x,y);
nx=1;
ny=25;
fx=sqrt(2/ax)*sin(pi*nx*x/ax);
fy=sqrt(2/ay)*sin(pi*ny*y/ay);
Z=zeros(p);
for i=1:p+1
    for j=1:p+1
        Z(i,j)=fx(i)*fy(j);    
    end
end
%Z=sqrt(4/(ax*ay))*sin(pi*nx*X./ax)*sin(pi*ny*Y./ay);

surf(X,Y,Z);
colormap(gray);
shading interp;
colorbar;
xlabel('x');
ylabel('y');
zlabel('Волновая функция');