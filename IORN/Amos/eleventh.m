function fact = eleventh (N)
    if N~=round(N) || N<0
        disp('error');
        fact = [];
    elseif N == 0
        fact = 1;
    else
        fact = N * eleventh (N-1);
    end
end