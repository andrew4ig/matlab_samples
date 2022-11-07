function crosspro = fourteenth (u,v)
    if length(u) == 2
        crosspro = [0 0 u(1)*v(2)-u(2)*v(1)];
    elseif length(u) == 3
        crosspro = [u(2)*v(3)-u(3)*v(2) u(3)*v(1)-u(1)*v(3) u(1)*v(2)-u(2)*v(1)];
    end
end