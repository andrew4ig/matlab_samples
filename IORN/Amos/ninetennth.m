function TriCirc = ninetennth (A,B,C)
    kAB=(B(2)-A(2))/(B(1)-A(1));
%     kAB= -1/kAB;
    kBC=(C(2)-B(2))/(C(1)-B(1));
%     kBC=-1/kBC;
    D=(A+B)/2;
    E=(B+C)/2;
    xo=((D(2)-E(2))*kBC*kAB+D(1)*kBC-E(1)*kAB)/(kBC-kAB);
    yo=D(2)-(xo-D(1))/kAB;
    [~,R]=sixteenth(A, B,C);
    plot([A(1), B(1), C(1), A(1)], [A(2), B(2), C(2), A(2)]);
    hold on;
    axis equal
    text( xo, yo, 'O');
    text( A(1), A(2), 'A');
    text( B(1), B(2), 'B');
    text( C(1), C(2), 'C');
    text( D(1), D(2), 'D');
    text( E(1), E(2), 'E');
    fimplicit(@(x,y) y-D(2)+(x-D(1))*kAB, [xo-R xo+R yo-R yo+R]);
    fimplicit(@(x,y) y-E(2)+(x-E(1))/kBC, [xo-R xo+R yo-R yo+R]);
    X=linspace(xo-R, xo+R);
    Y1=yo+(sqrt(R^2-(X-xo).^2));
    Y2=yo-(sqrt(R^2-(X-xo).^2));
    plot(X,Y1,X,Y2);
    hold off;
    grid on;
end