function ellipsePC = twentyth (xc, yc, a,b)
    fimplicit(@(x,y) ((x-xc)./a).^2+((y-yc)./b).^2-1', [xc-1.2*a xc+1.2*a yc-1.2*b yc+1.2*b]);
    ylabel('${y}$','Interpreter', 'latex');
    xlabel('${x}$','Interpreter', 'latex');
    axis equal;
    grid on;
    text(xc, yc,'$.(x_0,y_0)$','Interpreter', 'latex');
end