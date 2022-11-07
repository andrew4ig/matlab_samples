function z = myfun(A,B,q)
z=A*(1/B)*q-5;
end;
 
q=10
ezplot(@(x,y)myfun(x,y,q));
 