function T = polyadd (A, B, c)
    a=length(A); b=length(B);
    if(c=='add')
        c=1;
    elseif(c=='sub')
        c=-1;
    else
        fprintf('ERROR!');
    end
    
    if(a>b)
        T=zeros(1,a);
        T(1,a-b+1:a)=B*c;
        T=T+A;
    else
        T=zeros(1,b);
        T(1,b-a+1:a)=A;
        T=T+B*c;
    end
end