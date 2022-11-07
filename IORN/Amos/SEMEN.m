A1=@(x) 0.9^(32-x)*0.1^(x-20)/facttt(300-x,280)/facttt(20,20+x);

A2=0;
for i=20:32
   A2=A2+A1(i); 
end
A2=A2*facttt(269,300)/factorial(20)*0.1^20*0.9^268;

A2