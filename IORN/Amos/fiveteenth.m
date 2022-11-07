function S = fiveteenth (A,B,C)
    AB=B-A;
    AC=C-A;
    S=0.5*sqrt(sum(fourteenth(AB, AC).^2));
end