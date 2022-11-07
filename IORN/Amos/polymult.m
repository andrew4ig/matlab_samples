function T = polymult (A, B)
    a=length(A); b=length(B);
    n=(a-1)+(b-1)+1;
    T=zeros(1,n);
    for i=1:a
       for j=1:b
          k=n-(i-1)-(j-1);
          T(k)=T(k)+A(a-i+1)*B(b-j+1);           
       end
    end
end