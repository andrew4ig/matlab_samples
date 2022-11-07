function binary = eighteenth (d)
    if(d>2^16-1||d~=round(d)||d<0)
        fprintf('Error: number must be positive integer less than 2^16-1\n');
    else
        binary = 0; i=0;
        fprintf('%5i = ', d);
        while d~=0
            binary = binary + mod(d, 2)*10^i;
            d=floor(d/2);
            i=i+1;
        end
        fprintf('%17i\n', binary);
    end
end