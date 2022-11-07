function P = twsecond (M, N)
    if(M>=N)
        fprintf('\nN must be greater than M');
    elseif (M~=round(M) || N~=round(N) || M<0 || N<0)
        fprintf('\nN and M must be positive integers');
    else
       n=1;P=zeros(1);
       for i=M:N
          for j=2:round(i/2)+1
             if(mod(i,j)==0)
                 c=1;
                 break
             else
                c=0;
             end
          end
          if c==0
              n=n+1;
              P(n)=i;
          end
       end 
    end
end