function TWC = nineth(Ta,V)
    C1=35.74;
    C2=0.6215;
    C3=-35.75;
    C4=0.4275;
    TWC=round(C1+C2*Ta+C3*V^0.16+C4*Ta*V^0.16);
end