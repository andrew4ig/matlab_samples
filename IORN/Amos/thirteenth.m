function unitvec = thirteenth (A, B)
    AB=B-A;
    unitvec=AB/sqrt(sum(AB.^2));
end