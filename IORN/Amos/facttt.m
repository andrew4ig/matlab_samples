function N = facttt (a,b)
    N=1;
    if(a==b)
        N=1;
    else
        for i=a:b
            N=N*i;
        end
    end
end