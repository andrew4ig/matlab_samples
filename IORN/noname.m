A0 = 3; k = 5; w = 15;  dx = 0.1;
x = -10:0.1:10;
y = A0 * cos(k * x);
fig = figure;
a = axes(fig);
p = plot(a, x, y);
xlim([0 10]);
ylim([-A0*2 A0*2]);
time = 0; dt = 0.015;
% imitates time flow and calculates wave sum for each t-iteration
while true
    time = time + dt;
    p.XData = p.XData + dx;
    xlim([p.XData(1) p.XData(end)+1]);
    j=0;
    amplitude = zeros(1, length(p.XData));
    for X = x
        j=j+1;
        fun = @(k) exp(1i*(k * X - w * time));
        amplitude(j) = A0 * integral(fun, k-0.5,k+0.5);
    end
    p.YData = real(amplitude);
    drawnow
end
