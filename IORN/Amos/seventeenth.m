function circlePC = seventeenth (c, p)
    R=sqrt(sum((c-p).^2));
    fimplicit(@(x,y) (x-c(1)).^2+(y-c(2)).^2-R^2', [c(1)-1.2*R c(1)+1.2*R c(2)-1.2*R c(2)+1.2*R]);
    ylabel('${y}$','Interpreter', 'latex');
    xlabel('${x}$','Interpreter', 'latex');
    axis equal;
    grid on;
    text(c(1), c(2),'$.C$','Interpreter', 'latex');
    text(p(1), p(2),'$.P$','Interpreter', 'latex');
    hold on;
    plot([c(1) p(1)], [c(2), p(2)]);
    hold off;
end