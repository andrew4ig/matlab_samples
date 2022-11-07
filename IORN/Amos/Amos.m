clear
close all
clc
format SHORTG 
format compact
Chapter=input('Chapter №');
while(Chapter<=11&&Chapter>0)
    if Chapter==1
        fprintf('1.a) (22+5.1^2)/(50-6.3^2)=%f5\n', (22+5.1^2)/(50-6.3^2));
        fprintf('1.b) 44/7+8^2/5-99/3.9^2=%f5\n', 44/7+8^2/5-99/3.9^2);
        fprintf('2.a) sqrt(41^2-5.2^2)/(e^5-100.53)=%f\n', sqrt(41^2-5.2^2)/(exp(1)^5-100.53));
        fprintf('2.b) 132^1/3+ln(500)/8=%f\n', nthroot(132,3)+log(500)/8);
        fprintf('3.a) (14.8^3-6.3^2)/(sqrt(13)+5)^2=%g\n', (14.8^3-6.3^2)/(sqrt(13)+5)^2);
        fprintf('3.b) 45(288/9.3-4.6^2)-1065e^-1.5=%g\n', 45*(288/9.3-4.6^2)-1065*exp(1)^-1.5);
        fprintf('4.a) (24.5+64/3.5^2+8/3*12.5^2)/(sqrt(76.4)-28/15)=%g\n', (24.5+64/3.5^2+8/3*12.5^2)/(sqrt(76.4)-28/15))
        fprintf('4.b) (5.9^2-2.4^2)/3+(log_10(12890)/e^0.3)^2=%g\n', (5.9^2-2.4^2)/3+(log(12890)/log(10)/exp(1)^0.3)^2);
        fprintf('5.a) cos(7pi/9)+tan(7pi/15)sin(15°)=%g\n', cos(7*pi/9)+tan(7*pi/15)*sind(15));
        fprintf('5.b) sin^2(80°)-(cos14°sin80°)^2/0.18^1/3=%g\n', sind(80)^2-(cosd(14)*sind(80))^2/nthroot(0.18,3));
        x=6.7;
        fprintf('\n6.x=%g\na)0.01x^5-1.4x^3+80x+16.7=%g\nb)(x^3+e^x-51/x)^0.5=%g\n', [x,0.01*x^5-1.4*x^3+80*x+16.7, sqrt(x^3+exp(1)^x-51/x)]);
        clear x;
        t=3.9;
        fprintf('\n7.t=%g\na)56t-9.81t^2/2=%g\nb)14e^(-0.1t)sin(2pit)=%g\n', [t 56*t-9.81*t^2/2 14*exp(1)^(-0.1*t)*sin(2*pi*t)]);
        clear t;

        x=5.1; y=4.2;
        fprintf('\n8.x=%g, y=%g\na)3/4*xy-7x/y^2+(xy)^0.5=%g\nb)(xy)^2-(x+y)/(x-y)^2+((x+y)/(2x-y))^0.5=%g\n', [x y 3/4*x*y-7*x/y^2+sqrt(x*y) (x*y)^2-(x+y)/(x-y)^2+sqrt((x+y)/(2*x-y))])
        clear x y;

        a=12; b=5.6; c=3*a/b^2; d=(a-b)^2/c;
        fprintf('\n9.a=%g b=%g c= 3*a/b^2=%g d=(a-b)^2/c=%g\na)a/b+(d-c)/(d+c)-(d-b)^2=%g\nb)e^((d-c)/(a-2b))+ln(|c-d+b/a|)=%g\n',[a b c d a/b+(d-c)/(d+c)-(d-b)^2 exp(1)^((d-c)/(a-2*b))+log(abs(c-d+b/a))])
        clear a b c d

        syms x positive
        r=24;
        a1=solve(4*pi*r^3/3==x*x/2*x/4, 'Real',true);
        a2=solve(4*pi*r^2==2*(x*x/2+x*x/4+x/2*x/4), x, 'Real', true);
        fprintf('\n10. r=%g\na) a=%g b=a/2=%g c=a/4=%g\nb)a=%g b=a/2=%g c=a/4=%g\n',[r a1 a1/2 a1/4 a2 a2/2 a2/4]);
        clear x r a1 a2

        a=11; b=9;
        L=0.5*sqrt(b^2+16*a^2)+b^2/8/a*log((4*a+sqrt(b^2+16*a^2))/b);
        fprintf('\n11.L= 0.5*sqrt(b^2+16*a^2)+b^2/(8*a)*log((4*a+sqrt(b^2+16*a^2))/b)= 0.5*sqrt(%g^2+16*%g^2)+%g^2/(8*%g)*log((4*%g+sqrt(%g^2+16*%g^2))/%g)=%g\n', [b a b a a b a b L]);
        clear a b L

        x=pi/12;
        fprintf('\n12.x=pi/12\na) sin(5x) = 5sin(x)- 20sin^3(x)+16sin^5(x)\nsin(5x) = sin(5*%g) = sin(%g) = %g\n5sin(x)- 20sin^3(x)+16sin^5(x) = 5sin(%g)- 20sin^3(%g)+16sin^5(%g) = 5*%g-20*%g^3+16*%g^5 = %g\12b) sin^2(x)cos^2(x) = (1-cos(4x))/8\nsin^2(x)cos^2(x) = sin^2(%g)cos^2(%g) = (%g)^2*(%g)^2 = %g\n(1-cos(4x))/8 = (1-cos(4*%g))/8 = (1-cos(%g))/8=%g\nTrue.\n', [x 5*x sin(5*x) x x x sin(x) sin(x) sin(x) 5*sin(x)-20*sin(x)^3+16*sin(x)^5 x x sin(x) cos(x) sin(x)^2*cos(x)^2 x 4*x (1-cos(4*x))/8])
        clear x

        x=24;
        fprintf('\n13. x = 24°\na) tan(3x) = (3tan(x)-tan(x)^3)/(1-3tan(x)^2)\ntan(3x) = tan(3*%g) = tan(%g) = %g \n (3tan(x)-tan(x)^3)/(1-3tan(x)^2) = (3tan(%g)-tan(%g)^3)/(1-3tan(%g)^2) = (%g-%g)/(1-%g) = %g/%g = %g\n%g=%g\nTrue.\n', [x 3*x tand(3*x) x x x 3*tand(x) tand(x)^3 3*tand(x)^2 3*tand(x)-tand(x)^3 1-3*tand(x)^2 (3*tand(x)-tand(x)^3)/(1-3*tand(x)^2) tand(3*x) (3*tand(x)-tand(x)^3)/(1-3*tand(x)^2)]);
        clear x

        a=pi/6; b=3*pi/8;
        fprintf('\n14.a=pi/6 b=3pi/8\nsin(a)+sin(b)=2sin((a+b)/2)cos((a-b)/2)\nsin(a)+sin(b)= sin(%g)+sin(%g)= %g+%g= %g\n 2sin((a+b)/2)cos((a-b)/2)= 2sin((%g+%g)/2)cos((%g-%g)/2)= 2sin(%g)cos(%g)= 2*%g*%g= %g\n%g=%g\nTrue.', [a b sin(a) sin(b) sin(a)+sin(b) a b a b (a+b)/2 (a-b)/2 sin((a+b)/2) cos((a-b)/2) 2*sin((a+b)/2)*cos((a-b)/2) 2*sin((a+b)/2)*cos((a-b)/2) sin(a)+sin(b)]);
        clear a b

        x=linspace(pi/3,3*pi/2);
        y=x.*sin(0.6*x);
        fprintf('\n15./xsin(ax)dx=sin(ax)/a^2-xcos(ax)/a\n/xsin(0.6x)dx = %g\nsin(0.6x)/0.6-xcos(0.6x)/0.6 (3pi/2_pi/3) = (%g -%g)- (%g -%g)=%g\n',[sum(y)*(x(2)-x(1)) sin(0.6*3*pi/2)/0.6^2 3*pi/2*cos(0.6*3*pi/2)/0.6 sin(0.6*pi/3)/0.6^2 pi/3*cos(0.6*pi/3)/0.6 sin(0.6*3*pi/2)/0.6^2-3*pi/2*cos(0.6*3*pi/2)/0.6-sin(0.6*pi/3)/0.6^2+pi/3*cos(0.6*pi/3)/0.6]);
        clear x y

        a=5.3; g=42; b=6; c=sqrt(a^2+b^2-2*a*b*cosd(g)); alp=acosd((c^2+b^2-a^2)/(2*c*b)); betha=acosd((a^2+c^2-b^2)/(2*a*c));
        fprintf('\n16. a=5.3 g=42° b=6 \na)c^2=a^2+b^2-2abcos(g)\nc^2=%g^2+%g^2-2*%g*%g*cos(%g°)\nc^2=%g => c=%g', [a b a b g a^2+b^2-2*a*b*cosd(g) sqrt(a^2+b^2-2*a*b*cosd(g))])
        fprintf('\nb)cos(apl)=(c^2+b^2-a^2)/(2*c*b)\nalp=%g°\ncos(betha)=(a^2+c^2-c^2)/(2*a*c)\nbetha=%g°',[acosd((c^2+b^2-a^2)/(2*c*b)) acosd((a^2+c^2-c^2)/(2*a*c))])
        fprintf('\nc)alp+betha+g=%g+%g+%g=%g\n', [alp betha g alp+betha+g]);
        clear a b c alp betha g

        a=5; g=25; b=7; c=sqrt(a^2+b^2-2*a*b*cosd(g)); alp=acosd((c^2+b^2-a^2)/(2*c*b)); betha=acosd((a^2+c^2-b^2)/(2*a*c));
        fprintf('\n17. a=5 g=25° b=7\na)c^2=a^2+b^2-2abcos(g)\nc^2=%g^2+%g^2-2*%g*%g*cos(%g°)\nc^2=%g => c=%g', [a b a b g a^2+b^2-2*a*b*cosd(g) sqrt(a^2+b^2-2*a*b*cosd(g))])
        fprintf('\nb)cos(apl)=(c^2+b^2-a^2)/(2*c*b)\nalp=%g°\ncos(betha)=(a^2+c^2-c^2)/(2*a*c)\nbetha=%g°',[acosd((c^2+b^2-a^2)/(2*c*b)) acosd((a^2+c^2-b^2)/(2*a*c))])
        fprintf('\nc)(a-b)/(a+b)=tan((alp-betha)/2)/tan((alp+betha)/2)\n(a-b)/(a+b)=%g\ntan((alp-betha)/2)/tan((alp+betha)/2)=%g\n', [(a-b)/(a+b) tand((alp-betha)/2)/tand((alp+betha)/2)]);
        clear a b c alp betha g

        L=4; tetha=35;
        fprintf('\n18. L=4 tetha=35 \nV=V1+V2=pi*r^2*h/3+4/3*pi*r^3/2=\n=pi*(L*sind(tetha/2))^2*L*cosd(tetha/2)/3+2/3*pi*(L*sind(tetha/2))^3=\n=pi/3*(L*sind(tetha/2))^2*(L*cosd(tetha/2)+2*L*sind(tetha/2))=%g\n',pi/3*(L*sind(tetha/2))^2*(L*cosd(tetha/2)+2*L*sind(tetha/2)))
        clear L tetha

        a=48; g=83; b=34; 
        c=sqrt(a^2+b^2-2*a*b*cosd(g));
        s=(a+b+c)/2;
        r=a*b*c/(4*sqrt(s*(s-a)*(s-b)*(s-c)));
        fprintf('\n19. a=48 g=83° b=34 \na)c^2=a^2+b^2-2abcos(g)\nc^2=%g^2+%g^2-2*%g*%g*cos(%g°)\nc^2=%g => c=%g', [a b a b g a^2+b^2-2*a*b*cosd(g) sqrt(a^2+b^2-2*a*b*cosd(g))])
        fprintf('\nb)r=a*b*c/(4*sqrt(s*(s-a)*(s-b)*(s-c)))=%g\n',r)
        clear a b c s g r

        xa=2; ya=-3; za=1;
        x0=-4; y0=-2; z0=-3;
        a=0.6; b=0.5; c=0.7;
        d0=sqrt((xa-x0)^2+(ya-y0)^2+(za-z0)^2);
        d=d0*sin(acos(((xa-x0)*a+(ya-y0)*b+(za-z0)*c)/(d0*sqrt(a^2+b^2+c^2))));
        fprintf('\n20.d=d0*sin(acos(((xa-x0)*a+(ya-y0)*b+(za-z0)*c)/(d0*sqrt(a^2+b^2+c^2))))\nd0=sqrt((xa-x0)^2+(ya-y0)^2+(za-z0)^2)=%g\nd=d0*sin(acos(((xa-x0)*a+(ya-y0)*b+(za-z0)*c)/(d0*sqrt(a^2+b^2+c^2))))=%g\n',[d0 d]);
        clear xa ya za x0 y0 z0 a b c d0 d

        a=16; b=11;
        c=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)));
        fprintf('\n21.C=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)))=%g\n', c);
        clear a b c 

        Total=315; capacity=37;
        N=capacity*ceil(Total/capacity)-Total;
        fprintf('\n22.Amount of free seats: %g\n',N);
        clear Total capacity N

        Total=739; capacity=54;
        N=Total-capacity*fix(Total/capacity);
        fprintf('\n23.Amount of unpacked apples: %g\n',N);
        clear Total capacity N

        Number=316501.673;
        Num1=round(Number/100)*100;
        Num2=round(Number/1000)*1000;
        fprintf('\n24.Number=%g\nrounded to hundreds: %g\nrounded to thousands:%g\n', [Number Num1 Num2])
        clear Number Num1 Num2

        V0=14; R1=120.6; R2=119.3; R3=121.2; R4=118.8;
        V=V0*((R1*R3-R2*R4)/((R1+R2)*(R3+R4)));
        fprintf('\n25.V0=14; R1=120.6; R2=119.3; R3=121.2; R4=118.8;\nV=V0*((R1*R3-R2*R4)/((R1+R2)*(R3+R4)))=%g\n', V)
        clear V0 R1 R2 R3 R4 V

        L=0.15; R=14; C=2.6e-6;
        f=1/(2*pi)*sqrt(1/(L*C)-(R/L)^2);
        fprintf('\nL=0.15; R=14; C=2.6e-6;\n26.f=1/(2*pi)*sqrt(1/(L*C)-(R/L)^2)=%g\n', f)
        clear L R C f

        C=@(n,r)factorial(n)/factorial(r)/factorial(n-r);
        fprintf('\n27.a)C(49,6)=%g\nb)C(6,2)C(43,4)/C(49,4)=%g',[C(49,6),C(6,2)*C(43,4)/C(49,6)])
        clear C

        fprintf('\n28.a)log4(0.085)=log(0.085)/log(4)=%g\nlog6(1500)=log10(1500)/log10(6)=%g\n',[log(0.085)/log(4) log10(1500)/log10(6)])

        R1=120; R2=220; R3=75; R4=130;
        Req=1/(1/R1+1/R2+1/R3+1/R4);
        fprintf('\n29.R1=120; R2=220; R3=75; R4=130;\nReq=1/(1/R1+1/R2+1/R3+/1R4))=%g\n', Req);
        clear R1 R2 R3 R4 Req

        V0=36; R=2500; C=1600e-6; T=8;
        V=V0*(1-exp(-T/(R*C)));
        fprintf('\n30.V0=36; R=2500; C=1600e-6; T=8;\nV=V0*(1-exp(-T/(R*C)))=%g\n', V);
        clear V0 R C T V

        syms x t
        k=solve(exp(x)==2^(-1), 'Real', true);
        T=5730; ft=0.7745;
        time=solve(ft==exp(k*t/T), t, 'Real', true);
        fprintf('\n31.Time~%g(years)\n', round(time));
        clear x t k T ft time

        fprintf('\n32.a)GCD(91, 147)=%g\nb)GCD(555,962)=%g\n', [gcd(91, 147) gcd(555,962)]);

        syms x positive
        M1=9.5;
        M2=8.7;
        W1=solve(2/3*log(x)/log(10)-10.7==M1, 'Real',true);
        W2=solve(2/3*log(x)/log(10)-10.7==M2, 'Real',true);
        fprintf('\n33.M01/M02=%g\n', W1/W2);
        clear x M1 M2 W1 W2

        v=5000; c=3e8; L=2;
        d=L*(1-sqrt(1-v^2/c^2));
        fprintf('\n34.v=5000; c=3e8; L=2;\nd=L*(1-sqrt(1-v^2/c^2))=%g\n', d);
        clear v c L d

        n=5; P=80000; m1=365; m2=1; r=0.05;
        dif=P*(1+r/m1)^(n*m1)-P*(1+r/m2)^(n*m2);
        fprintf('\n35.n=5; P=80000; m1=365; m2=1;\ndif=%g\n', dif);
        clear n P m1 m2 dif r
        
        T1=79.5; T2=78; T0=98.6; Tout=69;
        Hours=21; Minutes=19;
        syms S1 S2 t k
        S1=T1==Tout+(T0-Tout)*exp(-k*t);
        S2=T2==Tout+(T1-Tout)*exp(-k);
        K=solve(S2,k);
        T=solve(subs(S1,k,K),t);
        T=floor(double(T)*60);
        Minutes=Minutes-rem(T,60);
        if Minutes<0
            Minutes=60+Minutes;
            Hours=Hours-1;
        end
        Hours=Hours-floor(T/60);
        fprintf('\n36.Time of the death %i:%i\n', [Hours Minutes]);
        clear S1 S2 T1 T2 T0 Tout t k T K Hours Minutes
        
        sg=12000; h=5; b=4; a=1.5;
        K=sg*sqrt(pi*a*h)*((1-a/(2*b)+0.326*(a/b)^2)/sqrt(1-a/b));
        fprintf('\n37.Koefficient of tension equals to %g\n', K);
        clear sg h b a K
        
        syms t
        N=20; N1=1e6;
        S1=N*exp(0.15*t);
        t1=double(solve(S1==2*N));
        t2=double(solve(S1==N1));
        fprintf('\n38.Time needed to double virus %.0g minutes\n', t1);
        fprintf('   Time needed to reach 1000000 %.0f minutes\n', t2);
        clear t N N1 t1 t2
        
        format rational
        fprintf('\n39.')
        disp(5/8+16/6);
        disp(1/3+11/13+2.7^2);
        format shortG
        
        nFac=@(n)sqrt(2*pi*n)*(n/exp(1))^n;
        d=@(n)(factorial(n)-nFac(n))/factorial(n);
        fprintf('\n40.20!=%g\n', nFac(20));
        fprintf('   Ponderousness=%g\n', d(20));
        clear nFac d
                
    elseif Chapter==2
        a=[8 10/4 12*1.4 51 tand(85) sqrt(26) 0.15];
        fprintf('1.%g %g %g %g %g %g %g\n', a);
        clear a

        a=[sqrt(15)*1e3 25/(10-6^2) log(35)/0.43 sind(65)/cosd(80)*129 cos(pi/20)^2];
        fprintf('2.%g %g %g %g %g \n',a);
        clear a

        a=[25.5; 14*tand(58)/(2.1^2+11); factorial(6); 2.7^4; 0.00552; pi/5];
        fprintf('3.%g %g %g %g %g %g\n', a);
        clear a

        a=[32/3.2^2; sind(35)^2; 6.1; log(29^2); 0.00552; log(29)^2; 133];
        fprintf('4.%g %g %g %g %g %g %g\n', a);
        clear a

        x=0.85; y=112.5;
        a=[y; y^x; log(y/x); x*y; x+y];
        fprintf('5.%g %g %g %g %g\n', a);
        clear x y a

        a=3.5; b=-6.4;
        x=[a a^2 a/b a*b sqrt(a)];
        fprintf('6.%g %g %g %g %g\n', x);
        clear a b x

        a=-1:6:-43;
        b=[8 10/4 12*1.4 51 tand(85) sqrt(26) 0.15];
        fprintf('7.%g %g %g %g %g %g %g %g\n%g %g %g %g %g %g %g\n', [a b]);
        clear a b

        s=linspace(96, 2, 11);
        fprintf('\n8.%g %g %g %g %g %g %g %g %g %g %g %g\n', s);
        clear s

        s=(26:-3.6:-10)';
        fprintf('\n9.%g %g %g %g %g %g %g %g %g %g %g %g\n', s);
        clear s

        s=linspace(-34,-7,9)';
        fprintf('\n10.%g %g %g %g %g %g %g %g %g %g\n', s);
        clear s

        Fives(1:5)=5;
        fprintf('\n11.%g %g %g %g %g', Fives);
        clear Fives

        Nines(1:9)=9;
        fprintf('\n12.%g %g %g %g %g %g %g %g %g', Nines);
        clear Nines

        a(6)=4.7;
        fprintf('\n13.%g %g %g %g %g %g', a);
        clear a

        b(6:8)=3.5;
        fprintf('\n14.%g %g %g %g %g %g %g %g', b);
        clear b

        b=[0 2 4 6 8 10 12 9 6 3 0];
        fprintf('\n15.%g %g %g %g %g %g %g %g %g %g %g', b);
        clear b

        a=2:3:17; b=3:4:15;
        c=[a b];
        fprintf('\n16.%g %g %g %g %g %g %g %g %g %g', c);
        clear a b c

        a=(2:3:17)'; b=(3:4:15)';
        c=[a; b];
        fprintf('\n17.%g %g %g %g %g %g %g %g %g %g', c);
        clear a b c

        vtA=8:7:71;
        vtB(1:4)=vtA(1:4);
        vtB(5:7)=vtA(8:10);
        fprintf('\n18.%g %g %g %g %g %g %g %g %g %g', vtB);
        clear vtA vtB

        vctC=5:4:49;
        vctC=reshape(vctC, [2, 6]);
        Codd=vctC(1,:);
        Ceven=vctC(2, :);
        fprintf('\n19.a)%g %g %g %g %g %g \n   b)%g %g %g %g %g %g', [Codd, Ceven]);
        clear vctC Codd Ceven

        vctD=0:3:27;
        vctDop=fliplr(vctD);
        fprintf('\n20. %g %g %g %g %g %g %g %g %g %g', vctDop);
        clear vctD vctDop

        A=[linspace(130, 10, 7); linspace(1, 12, 7); linspace(12, 72, 7)]; 
        fprintf('\n21.%g %g %g %g %g %g %g \n   %g %g %g %g %g %g %g\n   %g %g %g %g %g %g %g', [A(1,:),A(2,:),A(3,:)]);
        clear A

        B=[linspace(5, 5, 5); linspace(2, 2, 5); linspace(3, 3, 5)]';
        fprintf('\n22.%g %g %g\n   %g %g %g\n   %g %g %g\n   %g %g %g\n   %g %g %g', [B(1,:),B(2,:),B(3,:),B(4,:),B(5,:)]);
        clear B

        C=ones(2,5)*7;
        fprintf('\n23.%g %g %g %g %g\n   %g %g %g %g %g', [C(1,:),C(2,:)]);
        clear C

        D(1:3,5)=linspace(8,6,3);
        fprintf('\n24.%g %g %g %g %g\n   %g %g %g %g %g\n   %g %g %g %g %g', [D(1,:),D(2,:),D(3,:)]);
        clear D

        E(3:4,3:5)=reshape(linspace(5,0,6), 3, [])';
        fprintf('\n25.%g %g %g %g %g\n   %g %g %g %g %g\n   %g %g %g %g %g\n   %g %g %g %g %g', [E(1,:),E(2,:),E(3,:),E(4,:)]);
        clear E

        F(2:4,3:5)=[linspace(1,3,3);linspace(10,6,3);linspace(20,32,3)]';
        fprintf('\n26.%g %g %g %g %g\n   %g %g %g %g %g\n   %g %g %g %g %g\n   %g %g %g %g %g', [F(1,:),F(2,:),F(3,:),F(4,:)]);
        clear F

        a=[3 -1 5 11 -4 2];
        b=[7 -9 2 13 1 -2];
        c=[-2 4 -7 8 0 9];
        A=[a; b; c];
        B=A';
        fprintf('\n27.%g %g %g %g %g %g\n   %g %g %g %g %g %g\n   %g %g %g %g %g %g\n', [A(1,:),A(2,:),A(3,:)]);
        fprintf('\n   %2g %2g %2g\n   %2g %2g %2g\n   %2g %2g %2g\n   %2g %2g %2g\n   %2g %2g %2g\n   %2g %2g %2g', [B(1,:),B(2,:),B(3,:),B(4,:),B(5,:),B(6,:)]);
        clear a b c A B
        
        a=[3 -1 5 11 -4 2];
        b=[7 -9 2 13 1 -2];
        c=[-2 4 -7 8 0 9];
        A=[a(3:6); b(3:6); c(3:6)];
        B=[a(1:3)', b(1:3)', c(1:3)'];
        fprintf('\n28.a\n');
        disp(A);
        fprintf('\n28.b\n');
        disp(B);
        clear a b c A B
        
        a=[3 9 -0.5 3.6 1.5 -0.8 4];
        b=[12 -.8 6 2 5 3 -7.4];
        A=[a(3:6);a(4:7);b(2:5)];
        B=[a(2:7)',[b(1:3),b(5:7)]'];
        fprintf('\n29.a\n');
        disp(A);
        fprintf('\n29.b\n');
        disp(B);
        clear a b c A B

        A=[1 5 9 13 17];
        B=[1 5 9  1 5 9 13 17];
        C=[1 1; 5 5; 9 9;13 13; 17 17];
        D=[1 5 9 13 17;
           1 5 9 13 17];
        E=[1 5 9 13 17 1;
           1 5 9 13 17 5;
           1 5 9 13 17 9;
           1 5 9 13 17 13;
           1 5 9 13 17 17];
       fprintf('\n30.a\n');
       disp(A);
       fprintf('\n30.b\n');
       disp(B);
       fprintf('\n30.c\n');
       disp(C);
       fprintf('\n30.d\n');
       disp(D);
       fprintf('\n30.e\n');
       disp(E);
       clear A B C D E
      
       A=[-4 5 8 1 -0.2 -7];
       B=[6 -4 11 -4 5 8 1 -0.2 5 1];
       C=[19; 6; 8; 5];
       fprintf('\n31.a\n');
       disp(A);
       fprintf('\n31.b\n');
       disp(B);
       fprintf('\n31.c\n');
       disp(C);
       clear A B C
        
       A=[6 11 -4 -0.2 1 8; 5 6 5 8 1 11];
       B=[19 11 -4 5 6; 8 -4 5 11 -0.2; 5 -7 1 5 5]';
       fprintf('\n32.a\n');
       disp(A);
       fprintf('\n32.b\n');
       disp(B);
       clear A B
       
       A=[36:-2:26;
          24:-2:14; 
          12:-2:2];
       ha=A(1,:);
       hb=A(:,6);
       hc=[A(3,1:2), A(1,4:6)];
       fprintf('\n33.a) ha = \n');
       disp(ha);
       fprintf('\n33.b) hb = \n');
       disp(hb);
       fprintf('\n33.c) hc = \n');
       disp(hc);
       clear A ha hb hc
       
       A=1:18;
       B=reshape(A,3,[]);
       Ba=[B(:,1);B(:,3);B(:,5)];
       Bb=[B(2,2:5), B(3,:)];
       Bc=[B(1,3:5), B(3,2:4)];
       fprintf('\n34.a) Ba = \n');
       disp(Ba);
       fprintf('\n34.b) Bb = \n');
       disp(Bb);
       fprintf('\n34.c) Bc = \n');
       disp(Bc);
       clear A B Ba Bb Bc
      
       A=[1.5:0.5:5, 9.6:-0.5:6.1];
       D=reshape(A,[],4)';
       Da=[D(1,:),D(3,:)]';
       Db=[D(:,2);D(:,4)]';
       Dc=[D(1,1:2), D(2:4,2)', D(4,1:3)];
       fprintf('\n35.a) Da = \n');
       disp(Da);
       fprintf('\n35.b) Db = \n');
       disp(Db);
       fprintf('\n35.c) Dc = \n');
       disp(Dc);
       clear A D Da Db Dc
       
       E=[0, ones(1,5)*5; 0.1:0.2:0.7, 0.7:0.2:0.9; 12:-3:-3; 6:1:11];
       F=E(2:3,3:5);
       G=E(:,3:6);
       fprintf('\n36.a) F = \n');
       disp(F);
       fprintf('\n36.b) G = \n');
       disp(G);
       clear E F G
       
       H=[1.25:0.25:2.75;1:3 1:4;45:-5:15];
       G=[H(1,1:3), H(1,6:7); H(3,3:7)];
       K=[H(:,2)';H(:,3)';H(:,5)';H(:,7)'];
       fprintf('\n37.a) G = \n');
       disp(G);
       fprintf('\n37.b) K = \n');
       disp(K);
       clear H G K
       
       A=[1 13 16;
          3 15 18];
       B=[10 10 13 16;
          11 11 14 17;
          12 12 15 18];
       C=[1 4 7 10 13 16;
          2 5 8 11 14 17];
       D=[5 8;
          6 9];
       fprintf('\n38.a)\n');
       disp(A);
       fprintf('\n38.b)\n');
       disp(B);
       fprintf('\n38.c)\n');
       disp(C);
       fprintf('\n38.d)\n');
       disp(D);
       clear A B C D
      
       A=[2  12;
          10 20;
          18 32;
          29 44];
       B=[18 20 23 26 6 14 23 25 47];
       C=[0 0 0 0 0 0;
          0 0 0 0 0 0;
          0 0 0 0 32 44;
          0 0 0 0 35 47];
       fprintf('\n39.a)\n');
       disp(A);
       fprintf('\n39.b)\n');
       disp(B);
       fprintf('\n39.c)\n');
       disp(C);
       clear A B C
       
       v=[1 3 5 7 9 11 13 15 17 19 21 23];
       fprintf('\n 40. v = \n');
       disp(v);
       M=[1 7 13 19;
          3 9 15 21;
          5 11 17 23];
       fprintf('\n   M = \n');
       disp(M)
       M=[1 7 19;
          5 11 23];
       fprintf('\n   M = \n');
       disp(M)
       N=[1 1 1;
          1 1 1];
       fprintf('\n   N = \n');
       disp(N);   
       clear v M N
       
       A=[ones(2,2), zeros(2,2)];
       B=[eye(2), zeros(2,2), ones(2,2)];
       C=[ones(1,4);zeros(2,4)];
       fprintf('\n41.a)\n');
       disp(A);
       fprintf('\n41.b)\n');
       disp(B);
       fprintf('\n41.c)\n');
       disp(C);
       clear A B C
       
       A=[eye(2), zeros(2,2), ones(2,1)];
       B=[ones(2,4);[eye(2),zeros(2,2)]];
       C=[zeros(4,1), [ones(2,3);zeros(2,3)], [zeros(2,1); ones(2,1)]];
       fprintf('\n42.a)\n');
       disp(A);
       fprintf('\n42.b)\n');
       disp(B);
       fprintf('\n42.c)\n');
       disp(C);
       clear A B C
       
       A=eye(2);
       B=ones(2,2);
       C=zeros(2,2);
       D=[[A, B, C]; [C, B, A]]; 
       fprintf('\n43. D = \n');
       disp(D);
       clear A B C D
       
       A=ones(2,3);
       A=[[A' 0*A'];[0*A' A']]; 
       fprintf('\n44. A = \n');
       disp(A);
       clear A
       
    elseif Chapter==3

        x=[-3 -2 -1 0 1 2 3];
        y=x.^3-exp(0.5*x)+x;
        fprintf('1. %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g - x-values\n   %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g - y-values\n',[x y]);
        clear x y

        x=[1 2 3 4 5 6];
        y=(x+5).^3./(x.^2);
        fprintf('\n2. %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g - x-values\n   %5.3g %5.3g %5.3g %5.3g %5.3g %5.3g - y-values\n',[x y]);
        clear x y

        x=[1.5 2.5 3.5 4.5 5.5 6.5];
        y=(x+7).^4./(x+1)./sqrt(x);
        fprintf('\n3. %6g %6g %6g %6g %6g %6g - x-values\n   %6.5g %6.5g %6.5g %6.5g %6.5g %6.5g - y-values\n',[x y]);
        clear x y

        x=[20 30 40 50 60 70];
        y=(2*sind(x)+cosd(x).^2)./(sind(x).^2);
        fprintf('\n4. %3g° %3g° %3g° %3g° %3g° %3g° - x-values\n   %4.3g %4.3g %4.3g %4.3g %4.3g %4.3g - y-values\n',[x y]);
        clear x y

        s=[50 100 150 200 250 300];
        r=sqrt(s/pi)/2;
        V=4*pi*r.^3/3;
        fprintf('\n5. %4.3g %4.3g %4.3g %4.3g %4.3g %4.3g - s-values\n   %4.3g %4.3g %4.3g %4.3g %4.3g %4.3g - r-values\n   %4.3g %4.3g %4.3g %4.3g %4.3g %4.3g - V-values\n',[s r V]);
        clear s r V

        eps=8.85e-12;
        lmbd=1.7e-7;
        R=0.1;
        z=[0 2 4 6 8 10]*1e-2;
        E=lmbd/(2*eps)*(R.*z./((z.^2+R^2).^1.5));
        fprintf('\n6.a) %g %6.5g %6.5g %6.5g %6.5g %6.5g - z-values\n    %g  %6.5g %6.5g %6.5g %6.5g %6.5g - E-values\n',[z E]);
        z=(2:0.01:6)/100;
        E=lmbd/(2*eps)*(R.*z./((z.^2+R^2).^1.5));
        [Emax, dot]=max(E);
        zmax=z(dot);
        fprintf('\n6.b)%8g - Emax\n    %8g - zmax(m)\n', [Emax zmax]);
        clear eps lmbd R z E dot Emax zmax

        t=0:2:20;
        V=24; C=4000e-6; R=3800; 
        tao=R*C;
        i=V/R*exp(-t/tao);
        V=V*(1-exp(-t/tao));
        fprintf('\n7.')
        table(t', V', (i*1000)', 'variableNames', {'t(s)', 'V(V)', 'i(mA)'})
        clear t V C R tao i V T ans

        u=[23.5 -17 6];
        l1=(u(1)^2+u(2)^2+u(3)^2)^0.5;
        l2=sqrt(sum(u.^2));
        fprintf('\n8.a)l=%4.3g\n  b)l=%4.3g\n',[l1 l2])

        u1=[7 -4 -11];
        w=u1*l1/sqrt(sum(u1.^2));
        fprintf('\n9.w=%3.2gi%+3.2Gj%+3.2gk\n',w);
        clear l1 l2 u u1 w

        fprintf('\n10.a)[15/15 8/8 -6/-6]=[1 1 1]\n   b)[225 120 -90\n      120  64 -48\n      -90 -48  36]\n   c)15*15+8*8+(-6)*(-6)=325\n\n');

        u=[5 -6 9]; v=[11 7 -4];
        p1=sum(u.*v);
        p2=u*v';
        p3=dot(u,v);
        fprintf('\n11.a)sum(u.*v)=%g\n   b)u*v"=%g\n   c)dot(u,v)=%g\n',[p1 p2 p3]);
        clear u v p1 p2 p3

        v=[2 3 4 5 6];
        a=v*2;
        b=v.^3;
        c=v.^v;
        d=v./2;
        fprintf('\n12.a)a = [%g %3g %3g %3g %3g]\n   b)b = [%g %3g %3g %3g %3g]\n   c)c = [%g %3g %3g %3g %3g]\n   d)d = [%g %3g %3g %3g %3g]\n',[a b c d]);
        clear v a b c d 

        v=[8 6 4 2];
        a=v.^0;
        b=v-5;
        c=v.^(-0.5);
        d=v.^-2;
        fprintf('\n13.a)a = [%g %3g %3g %3g]\n   b)b = [%g %3g %3g %3g]\n   c)c = [%g %3g %3g %3g]\n   d)d = [%g %3g %3g %3g]\n',[a b c d]);
        clear v a b c d 

        x=[1 2 3 4 5]; y=x*2;
        z=(x+y).^2./(x-y);
        Z=x.*log(x.^2+y.^2)+sqrt(y.^3./(y-x).^2);
        fprintf('\n14.a)z = [%g %3g %3g %3g %3g]\n   b)z = [%.4g %6.4g %6.4g %6.4g %6.4g]\n',[z Z]);
        clear x y z Z

        r=1.6e3; s=14.2;
        t=[1 2 3 4 5];
        x=[0 2 4 6 8];
        y=[3 6 9 12 15];
        G=x.*t+r/s^2*(y.^2-x).*t;
        R=r*(-x.*t+y.*t.^2)/15-s.^2*(y-x.^2/2).*t;
        fprintf('\n15.a)G = [%6.5g %6.5g %6.5g %6.5g %6.5g]\n   b)R = [%6.5g %6.5g %6.5g %6.5g %6.5g]\n',[G R]);
        clear r s t x y G R

        o=[0 0 0];
        c=[-5 -2 11];
        b=[-7 9 6];
        a=[8 5 -4];
        roa=a-o; rob=b-o; roc=c-o;
        rab=-roa+rob;
        rac=-roa+roc;
        S=sqrt(sum(cross(rab,rac).^2))/2;
        fprintf('\n16.S = %6g\n',S);
        clear o c b a roa rob roc rab rac S

        o=[0 0 0];
        c=[-6 8 2];
        b=[1 3 6];
        a=[2 5 1];
        roa=a-o; rob=b-o; roc=c-o;
        rac=-roa+roc;
        V=dot(rob,cross(roa,rac));
        fprintf('\n17.V = %g\n',V);
        clear o c b a roa rob roc rab rac V

        u=[5 -2 4]; v=[-2 7 3]; w=[8 1 -3];
        Left=dot((u+v),cross((v+w),(w+u)));
        Right=2*dot(u,cross(v,w));
        fprintf('\n18.(u+v)*[(v+w)x(w+u)]=2u*(vxw)\n   (u+v)*[(v+w)x(w+u)]=%g\n   2u*(vxw)=%g \n   %g=%g\n   (u+v)*[(v+w)x(w+u)]=2u*(vxw)\n',[Left Right Left Right]);
        clear u v w Left Right

        r1=[6 -3 2];
        r2=[2 9 10];
        tetha=acos(dot(r1,r2)/(sqrt(sum(r1.^2))*sqrt(sum(r2.^2))));
        fprintf('\n19.tetha=%g\n', tetha);
        clear r1 r2 tetha

        R=14;
        B=[-R 0];
        C=[R 0];
        x=linspace(-R,R,102);
        x(1)=[];
        x(101)=[];
        y=sqrt(R^2-x.^2);
        A=[x' y'];
        rab=B-A;
        rac=C-A;
        alp=round(acosd(sum(rab.*rac,2)./(sqrt(sum(rab.^2,2)).*sqrt(sum(rac.^2,2)))));
        rab(:,3)=0;
        rac(:,3)=0;
        Alp=round(real(asind(sqrt(sum(cross(rab,rac,2).^2,2))./(sqrt(sum(rab.^2,2)).*sqrt(sum(rac.^2,2))))));
        fprintf('\n20.a)a=%g\n   b)a=%g', [mean(alp) mean(Alp)]);
        clear A alp Alp B C R rab rac x y

        g=9.81; v=162; a=70; t=1:5:31;
        x=v*cosd(a)*t;
        y=v*sind(a)*t-g*t.^2/2;
        r=sqrt(x.^2+y.^2);
        tet=atand(x./y);
        fprintf('\n21. %6.5g %6.5g %6.5g %6.5g %6.5g %6.5g %6.5g - r-values\n    %6.5g %6.5g %6.5g %6.5g %6.5g %6.5g %6.5g - tetha-values\n',[r tet]);
        clear g v a t x y r tet


        fprintf('\n22.a)n=5\n');
        format long
        n=0:5;
        a=2.^n./factorial(n);
        fprintf('S  = %15.15g\ne2 = %15.15g\n', [sum(a) exp(2)]);
        fprintf('\nb)n=10\n');
        n=0:10;
        a=2.^n./factorial(n);
        fprintf('S  = %15.15g\ne2 = %15.15g\n', [sum(a) exp(2)]);
        fprintf('\nc)n=50\n');
        n=0:50;
        a=2.^n./factorial(n);
        fprintf('S  = %15.15g\ne2 = %15.15g\n', [sum(a) exp(2)]);
        clear S n e2 a

        fprintf('\n23.a)n=10\n');
        n=1:10;
        a=(9/10).^n./n;
        fprintf('S  = %15.15g\nln10=%15.15g\n', [sum(a) log(10)]);
        fprintf('\nb)n=50\n');
        n=1:50;
        a=(9/10).^n./n;
        fprintf('S  = %15.15g\nln10=%15.15g\n', [sum(a) log(10)]);
        fprintf('\nc)n=100\n');
        n=1:100;
        a=(9/10).^n./n;
        fprintf('S  = %15.15g\nln10=%15.15g\n', [sum(a) log(10)]);
        clear S n ln10 a
        
        fprintf('\n24.a)n=5\n');
        n=1:5;
        a=2.^(-n);
        fprintf('S  = %15.15g\n', sum(a));
        fprintf('\nb)n=10\n');
        n=1:10;
        a=2.^(-n);
        fprintf('S  = %15.15g\n', sum(a));
        fprintf('\nc)n=40\n');
        n=1:40;
        a=2.^(-n);
        fprintf('S  = %15.15g\n', sum(a));
        clear S n ln10 a
        
        x=[1 0.5 0.1 0.01 0.001 0.00001 0.0000001];
        a=(cos(2*x)-1)./(cos(x)-1);
        fprintf('\n25.(cos(2*x)-1)/(cos(x)-1)\n');
        H=table(x', a', 'variableNames', {'x', 'lim'});
        format long
        disp(H);
        format shortG
        clear x a H
        
        x=[2 1.5 1.1 1.01 1.001 1.00001 1.0000001];
        a=(x.^(1/3)-1)./(x.^(1/4)-1);
        fprintf('\n26.(x^(1/3)-1)/(x^(1/4)-1)\n');
        H=table(x', a', 'variableNames', {'x', 'lim'});
        format long
        disp(H);
        format shortG
        clear x a H
    
        fprintf('\n27.Consumption of water\n');
        P=10:10:200;
        Q=1020*sqrt(P).*(1-0.01*sqrt(P));
        H=table(P', Q', 'variableNames', {'P', 'Q'});
        disp(H);
        clear P Q H
        
        R=0.08206; T=300; a=1.39; b=0.0391; n=1;
        V=0.1:0.02:1;
        P1=n*R*T./V;
        P2=n*R*T./(V-n*b)-n^2*a./V.^2;
        dp=(P1-P2)./P2*100;
        [~,n]=max(dp);
        fprintf('\n28.Volume of max inaccuracy %g\n', V(n));
        H=table(V', P1', P2', dp', 'variableNames', {'V', 'P_{ideal}', 'P_{Wqqls}', 'dP/P'});
        disp(H);
        clear R T a b n V P1 P2 dp H
        
        A=[1 -3 5; 2 2 4; -2 0 6];
        B=[0 -2 1; 5 1 -6; 2 7 -1];
        C=[-3 4 1; 0 8 2; -3 5 3];
        fprintf('\n29.Matrixes\n a)A + B = \n');
        disp(A+B);
        fprintf('\n  B + A = \n');
        disp(B+A);
        fprintf('\n b)(A + B) + C = \n');
        disp((A + B) + C);
        fprintf('\n  A + (B + C) = \n');
        disp( A + (B + C));
        fprintf('\n c)3*(A + B) = \n');
        disp(3*(A + B));
        fprintf('\n  3A+3B = \n');
        disp( 3*A+3*B);
        fprintf('\n d)A*B+A*C = \n');
        disp(A*B+A*C);
        fprintf('\n  A*(B + C) = \n');
        disp(A*(B + C));
        clear A B C
        
        fprintf('\n30.a)Когда матрицы перестановочны\n');
        fprintf('   b)Всегда\n');
        fprintf('   c)Когда произведение коммутативно\n');
        fprintf('   c)Всегда\n');
        
        A=randi([1 10], 4);
        fprintf('\n31.A*A = \n');
        disp(A*A);
        fprintf('\n  A.*A = \n');
        disp(A.*A);
        fprintf('\n  A/A = \n');
        disp(A\A);
        fprintf('\n  A./A = \n');
        disp(A./A);
        fprintf('\n  det(A) = \n');
        disp(det(A));
        fprintf('\n  inv(A) = \n');
        disp(inv(A));
        clear A
        
        A=magic(6);
        N=zeros(1,14);
        for i=1:6
           N(i)=sum(A(:,i)); 
        end
        for i=1:6
           N(i+6)=sum(A(i,:)); 
        end
        N(13)=sum(diag(A));
        N(14)=sum(diag(fliplr(A)));
        fprintf('\n32.magig square\n');
        disp(N);
        clear A N i
        
        A=[-4 3 1; 5 6 -2; 2 -5 4.5];
        B=[-18.2 -48.8 92.5]';
        a=det(A);
        B1=A;
        B1(:,1)=B;
        b1=(det(B1));
        B2=A;
        B2(:,2)=B;
        b2=(det(B2));
        B3=A;
        B3(:,3)=B;
        b3=(det(B3));
        x=b1/a; y=b2/a; z=b3/a;
        fprintf('\n33.Solution of system of equations\n x=%g y=%g z=%g', [x y z]);
        clear A B a B1 b1 B2 b2 B3 b3 x y z
        
        A=[ 2.5 -1 3 1.5 -2;
            3 4 -2 2.5 -1;
            -4 3 1 -6 2;
            2 3 1 -2.5 4;
            1 2 5 -3 4];
        B=[57.1 27.6 -81.2 -22.2 -12.2]';
        o=det(A);
        B1=A;
        B1(:,1)=B;
        b1=(det(B1));
        B2=A;
        B2(:,2)=B;
        b2=(det(B2));
        B3=A;
        B3(:,3)=B;
        b3=(det(B3));
        B4=A;
        B4(:,4)=B;
        b4=(det(B4));
        B5=A;
        B5(:,5)=B;
        b5=(det(B5));
        a=b1/o; b=b2/o; c=b3/o; d=b4/o; e=b5/o;
        fprintf('\n34.Solution of system of equations\n a=%g b=%g c=%g d=%g e=%g', [a b c d e]);
        clear a b c d e A B o b1 b2 b3 b4 b5 B1 B2 B3 B4 B5
        
        T1=[3 1 1 2 1];
        T2=[1 2 1 3 1];
        T3=[1 1 0 3 3];
        T4=[2 0 3 1 2];
        T5=[1 2 3 0 2];
        T=[128 118 112 112 104]*8;
        T1=min(floor(T./T1));
        T2=min(floor(T./T2));
        T3=min(floor(T./T3));
        T4=min(floor(T./T4));
        T5=min(floor(T./T5));
        n=[T1 T2 T3 T4 T5];
        fprintf('\n35. Quantity of packs:');
        disp(n);
        clear T1 T2 T3 T4 T5 T T0 n k i
        
        R1=4; R2=4; R3=6; R4=4; R5=3; R6=2; R7=2.5;
        V1=18; V2=18; V3=12; V4=28;
        R=[R1+R2+R3, -R2, -R3, 0;
           -R2, R3+R2+R6, 0, -R6;
           -R4, 0, R4+R5,    -R5;
           0, -R6, -R4, R7+R6+R4];
       V=-[-V1 V2 -V3 V4]';
       I=R\V;
       fprintf('\n36. Currents of counters:');
       disp(I');
       clear R1 R2 R3 R4 R5 R6 R7 V1 V2 V3 V4 R V I
        
       V1=40; V2=39; V3=36; 
       R1=16; R2=20; R3=10; R4=14; R5=8; R6=16; R7=10; R8=15; R9=6; R10=4;
       R=[R1+R2+R3, -R2, - R3, 0, 0;
          -R2, R2+R5+R4+R6,-R5,-R6, - R4;
          -R3, -R5,R3+R5+R7,-R7,0;
          0,-R6,-R7,R6+R7+R8+R9,-R8;
          0,-R4,0,-R8,R4+R8+R10];
       V=[-V1 0 -V2 V3 V1]';
       I=R\V;
       fprintf('\n37. Currents of counters:');
       disp(I');
       clear all
       
    elseif Chapter == 4

       T=90; R=90;
       HI=-42.379+2.04901523*T+10.14333127*R-0.22475541*R-6.83783e-3*T^2 - ...
            5.481717e-2*R^2+1.22874e-3*T^2*R+8.5282e-4*T*R^2-1.99e-6*T^2*R^2; 
       fprintf('\n1. Температура индекса тепла: %g\n', round(HI));
       clear T R HI

       F=1e5;  N=[5 6 7 8 9 10]; r=4.35/100;
       P=F*(r/12)./((1+r/12).^(12*N)-1);   
       fprintf('\n2.\n');
       disp([N; P]');
       clear F N r P

       k=log(2)*60/40; N0=1; t=(0:2:24);
       N=round(N0*exp(k*t));
       fprintf('\n3.\n');
       disp([t; N]');
       clear k N0 r P t N

       r2=[12 16 20 24 28]; r1=0.7*r2;
       S=pi^2*(r2.^2-r1.^2);
       V=0.25*pi^2*(r2.^2-r1.^2).*(r1+r2);
       fprintf('\n4.\n');
       disp([r2; r1; V; S]');
       clear r2 r1 S V

       W=500; L=120; h=50;
       x=10:20:110;
       T=W*L*sqrt(h^2+x.^2)/h./x;
       fprintf('\n5. L=%g, W=%g, h=%g\n',[L W h]);
       disp('           x             T  ')
       disp([x;T]');
       clear W L h x T

       %grades=input('ВВедите оценки:');
       grades=[92 74 53 61 100 42 80 66 71 78 91 85 79 68];
       lngth=length(grades);
       mn=mean(grades);
       sko=std(grades);
       med=median(grades);
       fprintf('\n6.Имеется %g оценок\n',lngth);
       fprintf('  Средняя оценка %g\n',mn);
       fprintf('  Стандартное отклонение %g\n',sko);
       fprintf('  Оценка посередине %g\n',med);
       clear grades med sko mn lngth

       h=4:4:40;
       tet=[2 2.9 3.5 4.1 4.5 5 5.4 5.7 6.1 6.4];
       R=h.*cosd(tet)./(1-cosd(tet));
       fprintf('\n7.          h        tetha            R\n');
       disp([h; tet; R]');
       clear h R tet

       T=13.3;
       syms x
       k=solve(exp(x*T)==2^(-1), 'Real', true);
       t=0:4:48;
       %A by A0 is just exp(kt)
       S=vpa(round(exp(k*t)*1000)/1000);
       fprintf('\n8.t      S\n');
       disp( [t; S]');
       clear k S t T x

       %L=input('Чему равна сумма ипотеки:');
       %N=input('На какое количество лет берется ипотека:');
       %r=input('Чему равна ваша процентная ставка:');
       r=4.5;
       L=250000;
       N=30;
       P=L*(r/1200)*(1+r/1200)^(12*N)/((1+r/1200)^(12*N)-1);
       fprintf('\n9.Ежемесячная оплата %2g летней ипотеки в сумме %6.6g с процентной ставкой %2.4g равна $%4.6g',[N L r P]);
       clear L N P r

       A=20000;   r=6.5;    P=391.32;
       n=6:6:60;
       B=A*(1+r/1200).^n-P*1200/r*((1+r/1200).^n-1);
       fprintf('\n10.         n            B          p\n');
       disp( [n; B; round((A-B)/A*1e4)/1e2]');
       clear A B r P n

       h=-500:500:1e4;
       p=29.921*(1-6.8753e-6*h);
       T=49.161*log(p)+44.932;
       fprintf('\n11.         h           T\n');
       disp( [h; T]');
       clear h p T

       s=600;
       a=10:0.1:120;
       h=2*s./a;
       alp=atan(h./(a/2));
       A=a+2./sin(alp)*2;
       H=h+2./cos(alp);
       S=A.*H/2;
       [O,n]=min(S);
       fprintf('\n12.При a=%g и h=%g Площадь будет равна %g и будет наименьшей\n', [a(n) h(n) O]);
       clear s a h S A H alp O n

       R=55; 
       a=8:0.25:55;
       b=sqrt(R^2-a.^2);
       A=a-8;
       B=b-20;
       S=A.*B;
       [O, n]=max(S);
       fprintf('\n13.При a=%g и и=%g Площадь будет равна %g и будет наибольшей\n', [a(n) b(n) O]);
       clear a b R A B S O n

       vrun=3; vswim=1; L=48; ds=30; dw=42;
       y=20:48;
       l1=sqrt(ds^2+y.^2);
       l2=sqrt(dw^2+(L-y).^2);
       t1=l1./vrun;
       t2=l2./vswim;
       t=t1+t2;
       [T, n]=min(t);
       fprintf('\n14.При y=%g время преодоления пути будет кратчайшим и будет равно T=%g\n',[y(n) T]);
       fprintf('sin(pfi)/sin(alp)=%g/%g=%g\n',[y(n)/l1(n) (L-y(n))/l2(n) y(n)/l1(n)/((L-y(n))/l2(n))]);
       fprintf('vrun/vswim=%g/%g=%g\n',[vrun vswim vrun/vswim]);
       clear vswim vrun L ds dw l1 l2 t1 t2 t T n y

       h=900; H=70; 
       x=50:0.5:1500;
       tet=atand(h./x)-atand((h-H)./x);
       [Tet, n]=max(tet);
       fprintf('\n15.При x=%g угол обзора цели будет равен %2.2g и будет максимальным\n',[x(n) Tet]);
       clear h H x tet Tet n

       M=20; b=0.25; t=0.01; a=0.25;
       alp=a/b; beta=(pi*alp)/2;
       C=sqrt(tan(beta)/beta)*(0.923+0.199*(1-sin(beta))^2)/cos(beta);
       sg=6*M/(t*b^2);
       K=C*sg*sqrt(pi*a);
       fprintf('\n16.Коэффициент интенсивности напряжения для балки шириной %gм и толщиной %gм\n с граничной трещиной %gм и приложенным моментом %g Н-м, есть %2.2g Па·(м)^0.5\n', [b t a M K]);
       clear a alp b beta C K M sg t

       v=50; ro=2000; h=500; Alp=90;
       syms x
       T=double(solve(Alp==v*x/ro, 'Real', true));
       t=linspace(0, T, 15);
       alp=v*t/ro;
       r=sqrt((ro*sin(alp)).^2+(h+ro-ro*cos(alp)).^2);
       tet=90-atand((ro*sin(alp))./(h+ro-ro*cos(alp)));
       fprintf('\n17.         t        tetha            r\n');
       disp( [t; tet; r]');
       clear v ro h Alp x T t alp r tet

       C=13.83; E=0.67;
       load ('constants.mat','kB');
       t=xlsread('4.18.xlsx')';
       sg=exp(C-E./(2*kB*t));
       fprintf('\n18.         t        sigma\n');
       disp( [t;sg]');
    
    elseif Chapter==5
       x=linspace(-1,5);
       y=(x.^2-3*x+7)./sqrt(2*x+5);
       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       plot(x,y);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('1','Interpreter', 'latex');
       grid on

       x=linspace(-4,9);
       y=(3*cos(x)-sin(x)).*exp(-0.2*x);
       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       plot(x,y);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('2','Interpreter', 'latex');
       grid on

       x=linspace(-4,4);
       y=(x.^2)./(2+sin(x)+x.^4);
       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       plot(x,y);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('3','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       syms f(x)
       f(x)=x^3-2*x^2-10*sin(x)^2-exp(0.9*x);
       df=diff(f,x);
       fplot(f(x),[-2 4]);
       hold on
       fplot(df,'--',[-2 4]);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('4','Interpreter', 'latex');
       legend('$f(x)$','${df(x)\over dx}$','Interpreter', 'latex');
       grid on

       x=linspace(-4,3);
       y=-3*x.^4+10*x.^2-3;
       figure('Units', 'normalized', 'OuterPosition', [1/12 1/3 1/3 0.5]);
       plot(x,y);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('5.1','Interpreter', 'latex');
       grid on
       x=linspace(-4,4);
       y=-3*x.^4+10*x.^2-3;
       figure('Units', 'normalized', 'OuterPosition', [7/12 1/3 1/3 0.5]);
       plot(x,y);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('5.2','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       fplot(@(x)(sin(3.*x)+cos(5.*x).^2).*exp(-0.2.*x), [-6 6]);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('6','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       syms f(x)
       f(x)=sin(x)^2*cos(2*x);
       df=diff(f,x);
       fplot(f(x),[-2 4]);
       hold on
       fplot(df,'--',[-2 4]);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('7','Interpreter', 'latex');
       legend('$f(x)$','${df(x)\over dx}$','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       fimplicit(@(x,y) (y-2.7).^2+(x-4.2).^2-7.5^2, [-5 15 -6 12])
       axis equal
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('8','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       t=linspace(-pi,pi);
       x=sin(t).*cos(t);
       y=1.5*cos(t);
       plot(x,y);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('9','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       t=linspace(0,2);
       x=cos(t).^3;
       y=sin(x).^3;
       u=sin(t);
       v=cos(t);
       plot(x,y,u,v);
       xlim([-2 2]);
       ylim([-2 2]);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('10','Interpreter', 'latex');
       legend('y(x)', 'v(u)','Interpreter', 'latex');
       grid on
       clear u v x y df f t

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       x1=linspace(-1,2.9,50);
       x2=linspace(3.1,7,50);
       y1=(x1.^2-5*x1-12)./(x1.^2-x1-6);
       y2=(x2.^2-5*x2-12)./(x2.^2-x2-6);
       plot(x1,y1,'k',x2,y2,'k');
       ylim([-20 20]);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('11','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       x1=linspace(-4,-2.01,50);
       x2=linspace(-1.99,4.99,50);
       x3=linspace(5.01,9,50);
       y1=(x1.^2+3*x1-5)./(x1.^2-3*x1-10);
       y2=(x2.^2+3*x2-5)./(x2.^2-3*x2-10);
       y3=(x3.^2+3*x3-5)./(x3.^2-3*x3-10);
       plot(x1,y1,'k',x2,y2,'k',x3,y3,'k');
       ylim([-20 20]);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('12','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       x=@(t) 3*t./(1+t.^3);
       y=@(t) 3*t.^2./(1+t.^3);
       t1=linspace(-30,-1.6, 500);
       t2=linspace(-0.6,40, 500);
       x1=x(t1); x2=x(t2);
       y1=y(t1); y2=y(t2);
       plot(x1,y1,'k',x2,y2,'b');
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('13','Interpreter', 'latex');
       grid on
       clear x y t1 t2 x1 x2 x3 y1 y2 y3

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       t=linspace(0,4*pi,500);
       x=13*cos(t)-2*cos(6.5*t);
       y=13*sin(t)-2*sin(6.5*t);
       plot(x,y);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('14','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       t=linspace(-4,3);
       x=(3.3-0.4*t.^2).*sin(t);
       y=(2.5-0.3*t.^2).*cos(t);
       plot(x,y);
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('15','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       tet=linspace(0,2*pi);
       r=2*sin(3*tet).*sin(tet);
       polarplot(tet, r);
       title('16','Interpreter', 'latex');
       grid on
       clear tet r

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       fimplicit(@(x,y) (y-3).^2/4^2+(x-2).^2/10^2-1,[-10 15 -4 10])
       axis equal
       ylabel('$y$','Interpreter', 'latex');
       xlabel('$x$','Interpreter', 'latex');
       title('17','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       X=[1850 1910 1950 1980 2000 2010];
       Y=[1.3 1.75 3 4.4 6 6.8];
       t=1850:2010;
       P=11.55./(1+18.7*exp(-0.0193*(t-1850)));
       plot(X,Y,'.k', t, P, 'b');
       xlim([1800 2200]);
       grid on;
       legend('real','analytic');
       ylabel('$P,mlrd$','Interpreter', 'latex');
       xlabel('$t, years$','Interpreter', 'latex');
       title('18','Interpreter', 'latex');
       clear X Y t P

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       eps=0.885e-12;
       Q=9.4e-16;
       q=2.4e-5;
       R=0.1;
       z=linspace(0,0.3);
       F=Q*q*z/(2*eps).*(1-z./sqrt(z.^2+R^2));
       [Fmax, zmax]=max(F);
       zmax=z(zmax);
       plot(z,F,'b',zmax,Fmax,'.', 'markersize',12);
       grid on;
       ylabel('$F,Newtons$','Interpreter', 'latex');
       xlabel('$z, meters$','Interpreter', 'latex');
       title('19','Interpreter', 'latex');
       text(zmax-0.04/2,Fmax-2e-11,'Max value')
       clear Q q R z F Fmax zmax eps

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       t=linspace(0,20);
       r=25+30*(1-exp(0.07*t));
       tet=2*pi*(1-exp(-0.2*t));
       polarplot(tet,r);
       title('20','Interpreter', 'latex');
       grid on

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       dr=diff(r);
       dr(100)=dr(99);
       dtet=diff(tet);
       dtet(100)=dtet(99);
       dt=t(5)-t(4);
       vr=dr/dt;
       vtet=r.*dtet/dt;
       v=sqrt(vr.^2+vtet.^2);
       ylabel('$V,m/s$','Interpreter', 'latex');
       xlabel('$t,s$','Interpreter', 'latex');
       plot(t,v);
       title('21','Interpreter', 'latex');
       grid on
       clear dr dtet dt vr vtet v r tet

       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       t=linspace(0, 5, 500);
       x=52*t-9*t.^2;
       y=125-5*t.^2;
       dx=diff(x);
       dy=diff(y);
       dt=t(9)-t(8);
       vx=dx/dt;
       vy=dy/dt;
       v=sqrt(vx.^2+vy.^2);
       [Vmin, n]=min(v);
       t(500)=[];
       plot(x,y, x(n), y(n),'.','MarkerSize', 12);
       hold on
       plot(t*16,3*v-70,16*t(n),3*Vmin-70,'.','MarkerSize', 12);
       ylabel('${V,m/s \over y,m}$','Interpreter', 'latex');
       xlabel('${t,s \over x,m}$','Interpreter', 'latex');
       title('22','Interpreter', 'latex');
       grid on
       clear t x y dx dy dt vx vy v Vmin n

       P=0:200;
       Q=1020*sqrt(P).*(1-0.01*sqrt(P));
       figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
       plot(P, Q);
       ylabel('$Q, g/min$','Interpreter', 'latex');
       xlabel('$P, pop$','Interpreter', 'latex');
       title('23','Interpreter', 'latex');
       grid on
       clear P Q

       figure('Units', 'normalized', 'OuterPosition', [1/4 0 1/2 1]);
       t=linspace(0,20,200);
       syms x(t)
       x(t)=(-3+4*t)*exp(-0.4*t);
       v=diff(x,t);
       a=diff(v,t);
        subplot(3,1,1)
       fplot(x(t),[0 20]);
       ylabel('$x, m$','Interpreter', 'latex');
       xlabel('$t, s$','Interpreter', 'latex');
       title('24.1 x=x(t)','Interpreter', 'latex');
       grid on
        subplot(3,1,2)
       fplot(v,[0 20]);
       ylabel('$v, m/s$','Interpreter', 'latex');
       xlabel('$t, s$','Interpreter', 'latex');
       title('24.2 v=v(t)','Interpreter', 'latex');
       grid on
       subplot(3,1,3)
       fplot(a,[0 20]);
       ylabel('$a, m/s^2$','Interpreter', 'latex');
       xlabel('$t, s$','Interpreter', 'latex');
       title('24.3 a=a(t)','Interpreter', 'latex');
       grid on
       clear x v a t
    
    elseif Chapter==6
       fprintf('\n1.a)12-4<5*3 - %g',12-4<5*3);
       fprintf('\n  b)8/4>6*3-4^2-3 - %g',8/4>6*3-4^2-3);
       fprintf('\n  c)-3<(8-12)+2*(5>18/6-4)^2 - %g',-3<(8-12)+2*(5>18/6-4)^2);
       fprintf('\n  d)(~5+~0)*6==3+3*~0 - %g\n',(~5+~0)*6==3+3*~0);

       a=-2; b=3; c=5;
       fprintf('\n2.a)a-b>a-c<b - %g',a-b>a-c&&a-c<b);
       fprintf('\n  b)-4<a<0 - %g',-4<a&&a<0);
       fprintf('\n  c)a-c<=b>a+c - %g',a-c<=b&&b>a+c);
       fprintf('\n  d)3*(c+a~=a/b-b)==(a+c)~=b - %g\n',3*(c+a~=a/b-b)==(a+c)&&(a+c)~=b);
       clear a b c

       v=[4 -1 2 3 1 -2 5 0];
       u=[5 -1 0 3 -3 2 1 5];
       fprintf('\n3.a)~~u - [%g %g %g %g %g %g %g %g]',~~u);
       fprintf('\n  b)v==~u - [%g %g %g %g %g %g %g %g]',v==~u);
       fprintf('\n  c)u==abs(v) - [%g %g %g %g %g %g %g %g]',u==abs(v));
       fprintf('\n  d)v>=u+v - [%g %g %g %g %g %g %g %g]\n',v>=u+v);

       fprintf('\n4.w=')
       w=(u<=v).*u;
       disp(w);
       clear u v w

       fprintf('\n5.a)-3&3 - %g',-3&3);
       fprintf('\n  b)-5<4&~0>-3 - %g',-5<4&~0>-3);
       fprintf('\n  c)-2&2>3|8/3 - %g',-2&2>3|8/3);
       fprintf('\n  d)-3<-1~0|5<4<3 - %g\n',-3<-1*(~0)||5<4&&4<3);

       M(3,5)=0;
       for i=1:3
           for j=1:5
                M(i,j)=i^j/(i+j);           
           end
       end
       fprintf('\n6.\n')
       disp(M);

       M(1,1:7)=1;
       M(1:7,1)=1;
       for i=2:7
           for j=2:7
                M(i,j)=M(i-1,j)+M(i,j-1);           
           end
       end
       fprintf('\n7.\n')
       disp(M);
       clear M

       BOS=[2.67 1.00 1.21 3.09 3.43 4.71 3.88 3.08 1.10 2.62 1.01 5.93];
       SEA=[6.83 3.63 7.20 2.68 2.05 2.96 1.04 0.00 0.03 6.71 8.28 6.85];
       BOSsum=sum(BOS);
       BOSmean=mean(BOS);
       SEAsum=sum(SEA);
       SEAmean=mean(SEA);
       BOScount=sum(double(BOS>BOSmean));
       SEAcount=sum(double(SEA>SEAmean));
       n=BOS<SEA;
       fprintf('\n8.a)BOSsum = %g, SEAsum = %g',[BOSsum SEAsum]);
       fprintf('\n    BOSmean = %g, SEAmean = %g',[BOSmean SEAmean]);
       fprintf('\n  b)BOScount = %g',BOScount);
       fprintf('\n    SEAcount = %g',SEAcount);
       fprintf('\n  c)n = %g',       sum(double(n)));
       fprintf('\n    mounths = [%g %g %g %g %g %g %g %g %g %g %g %g]',n);
       clear BOS BOScount BOSmean BOSsum SEA SEAcount SEAmean SEAsum n

       i=1;
       while(true)
           if(mod(i,13)==0&&mod(i,16)==0&&sqrt(i)>120)
               break
           end
           i=i+1;
       end
       fprintf('\n9.Такое число есть: %g',i);

       a=0; b=1; i=2;
       fprintf('\n10.Fibonaci sequence: %g %g',[a b]);
       while (i<20)
          c=a+b;
          a=b;
          b=c;
          fprintf(' %g' ,c);
          i=i+1;
       end
       clear a b c i j

       a=1; b=1; i=2;
       F(1)=a; F(2)=b;
       while (i<100)
          i=i+1;
          c=a+b;
          a=b;
          b=c;
          F(i)=c;
       end
       S1=sum(1./F(1:10));
       S2=sum(1./F(1:50));
       S3=sum(1./F(1:100));
       fprintf('\n11. n=10  Psi = %.15g', S1);
       fprintf('\n    n=50  Psi = %.15g', S2);
       fprintf('\n    n=100 Psi = %.15g', S3);
       clear a b c F i S1 S2 S3


    %    a=input('a = ');
    %    b=input('b = ');
    %    c=input('c = ');
       fprintf('\n12.a)');
       a=3;b=6;c=3;
       quadroots(a,b,c);
       fprintf('\n   b)');
       a=-3;b=4;c=-6;
       quadroots(a,b,c);
       fprintf('\n   c)');
       a=-3;b=7;c=5;
       quadroots(a,b,c);
       clear a b c

       format long
       fprintf('\n13. Pi = ');
       n=1:1e6;
       n=n.^(-2);
       pie=sqrt(6*sum(n));
       disp(pie);
       fprintf('    Pi = ');
       disp(pi);
       clear n pie

       fprintf('\n14. Pi = ');
       b=0; c=1;
       for i=1:40
           b=sqrt(2+b);
           c=c*b;
       end
       c=c/2^40;
       pie=2/c;
       disp(pie);
       fprintf('    Pi = ');
       disp(pi);
       clear pie c b ans

       A=randi([-10 10],1,20);
       B=abs(A);
       S=(A+B)/2;
       S=sum(S);
       fprintf('\n15.');
       disp(A);
       fprintf('Sum of positive elements: %g\n', S);
       clear A B S i


       A=randi([-10 10],1,20);
       fprintf('\n16. A=\n');
       disp(A);
       n=0;
       B=abs(A);
       while(sum(A-B)~=0)
           n=n+1;
           for i=1:20
                if(A(i)<0)
                    A(i)=randi([-10 10]);
                end
           end
           B=abs(A);
       end
       fprintf('     A=\n');
       disp(A);
       fprintf('Number of iterations: %g\n', n);
       clear A B n i

       n=randi([5 30]);
       A=randi([-10 10],1,n);
       clear n
       n=length(A);
       pos=0;neg3=0;
       for i=1:n
          if(A(i)>0)
              pos=pos+1;
          end
          if(mod(A(i),3)==0&&A(i)<0)
             neg3=neg3+1; 
          end
       end
       fprintf('\n17.Вектор имеет %g элементов. %g элементов положительны, и %g элементов отрицательны и делятся на 3\n', ...
           [n pos neg3]);
       disp(A);
       clear A n pos neg3 i

       format shortG
       x=[4.5 5 -16.12 21.8 10.1 10 -16.11 5 14 -3 3 2];
       for i=1:length(x)
           for j=1:(length(x)-i)
                if(x(j+1)<x(j))
                   temp=x(j);
                   x(j)=x(j+1);
                   x(j+1)=temp;
                   clear temp
                end
           end
       end
       fprintf('\n18.Sorted array:\n');
       disp(x);
       clear x i j

       n=0; M(3)=0;
       for a=1:50
           for b=1:50
               for c=1:50
                   if(a*a+b*b==c*c)
                       n=n+1;
                       M(n,1:3)=[a b c];
                   end
               end
           end
       end
       fprintf('\n19.Pethagorian numbers\n');
       disp(M);
       clear n a b c M

       n=1;P=zeros(1,500); P(1)=3;
       for i=2:500
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

       n=0; M(2)=0;
       for i=2:length(P)
         if(P(i)-P(i-1)==2)
             n=n+1;
             M(n,1:2)=[P(i) P(i-1)];
         end
       end
       fprintf('\n20.Paired primes\n');
%        for j=1:24
%            for i=1:500
%               if(M(24-j+1,500-i+1)==0)
%                 M(24+1-j,500-i+1)=[];
%               end
%            end
%        end
       disp(M);
       clear n c i j M

       n=0; M(500)=0;
       for i=3:25
         if(P(i-1)~=P(i)-2&&P(i+1)~=P(i)+2)
             n=n+1;
             M(n)=P(i);
         end
       end
       fprintf('\n21.Isolatedd primes\n');
       for i=1:length(M)
          if(M(500-i+1)==0)
            M(500-i+1)=[];
          end
       end
       disp(M);
       clear i M n P

       S=[31 70 92 5 47 88 81 73 51 76 80 90 55 23 43 98 36 87 22 61 19 69 26 82 89 99 71 59 49 64];
       a20=0; a40=0; a60=0; a80=0; a100=0;
       for i=1:length(S)
           if S(i)>0 &&S(i)<20
               a20=a20+1;
           elseif S(i)>19 &&S(i)<40
               a40=a40+1;
           elseif S(i)>39 &&S(i)<60
               a60=a60+1;
           elseif S(i)>59 &&S(i)<80
               a80=a80+1;
           elseif S(i)>79 &&S(i)<100
               a100=a100+1;
           end
       end
       fprintf('\n22.Оценки от 0 до 19  %10g студента\n', a20);        
       fprintf('   Оценки от 20 до 39 %10g студента\n', a40);    
       fprintf('   Оценки от 40 до 59 %10g студента\n', a60);
       fprintf('   Оценки от 60 до 79 %10g студента\n', a80);
       fprintf('   Оценки от 80 до 99 %10g студента\n', a100);
       clear S a20 a40 a60 a80 a100 i
   
    elseif Chapter == 7
        fprintf('\n1. y=(-0.2*x^2+7*x^2)*exp(-0.17*x)');
        fprintf('\na) y(-1.5)=%g, y(5)=%g\n', [first(-1.5) first(5)]);
        x=linspace(-2,6);
        y=first(x);
        figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
        plot(x,y);
        ylabel('$y$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('1.b','Interpreter', 'latex');
        grid on
        clear x y

        fprintf('\n2. r=4cos(4sin(tet))');
        fprintf('\na) yr(pi/6)=%g, r(pi/6)=%g\n', [second(pi/6) second(5*pi/6)]);
        x=linspace(0,pi*2);
        y=second(x);
        figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
        polarplot(x,y);
        title('2.b','Interpreter', 'latex');
        grid on
        clear x y

        fprintf('\n3. gml=235.2lkm');
        fprintf('\na) 5gml=%gkm', third(5));
        fprintf('\nb) 5.8gml=%gkm\n', third(5.8));

        fprintf('\n4. ft/dm3=27000kg/m3');
        fprintf('\na) 7860kg/m3=%gft/d3', fourth(7860));
        fprintf('\nb) 4730kg/m3=%gft/d3\n', fourth(4730));

        fprintf('\n5. 1 khp=1.69 fps');
        fprintf('\n   400 khp=%3g fps\n', fifth(400));

        fprintf('\n6. BSA=0.0077184W^0.425*H^0.75;');
        fprintf('\na) for m = 95, h = 1.87: BSA = %5g', sixth(95, 1.87));
        fprintf('\nb) for m = 61, h = 1.58: BSA = %5g\n', sixth(61, 1.58));

        r=20;
        figure('Units', 'normalized', 'OuterPosition', [1/3 1/3 1/3 0.5]);
        plot(linspace(0,2*r), seventh(r, 1.5*r, 2*r, linspace(0,2*r)));
        ylabel('$V$','Interpreter', 'latex');
        xlabel('$H$','Interpreter', 'latex');
        title('Volume of fuel bank','Interpreter', 'latex');
        grid on
        clear r

        W=@(r,d,t,y) pi^2*(2*r+d)*d*t*y;
        y=0.696; r=0.35; d=0.12; t=0.002;
        fprintf('\n8.Вес покрытия равен W=%g\n', W(r,d,t,y));
        clear W y r d t

        fprintf('\n9. TWC=C1+C2*Ta+C3*V^0.16+C4*Ta*V^0.16;;');
        fprintf('\na) for T = 35, V = 26: Twc = %g', nineth(35, 26));
        fprintf('\nb) for T = 10, V = 50: Twc = %g\n', nineth(10, 50));

        g=[3.7 3 3.3 2 0 4 1.3 4];
        h=[4 3 3 2 3 4 3 3];
        fprintf('\n10. GPA = %g\n', tenth(g,h));
        clear g h

        fprintf('\n11. n!=n*(n-1)*...*2*1')
        fprintf('\n  9!   =');disp(eleventh(9));
        fprintf('  8.5! =');disp(eleventh(8.5));
        fprintf('  0!   =');disp(eleventh(0));
        fprintf('  -5!  =');disp(-eleventh(5));

        fprintf('\n12. th=acos(BA*BC/(2|BA|*|BC|))');
        fprintf('\na) for A = (-5, 1, 6), B = (2.5, 1.5, -3.5), C = (2.3, 8, 1): th = %.2f°', twelwth([-5, 1, 6], [2.5, 1.5, -3.5], [2.3, 8, 1]));
        fprintf('\nb) for A = (-5.5, 0), B = (3.5, -6.5), C = (0, 7): th = %.2f°\n', twelwth([-5.5, 0], [3.5, -6.5], [0, 7]));

        fprintf('\n13. unitvec=AB/|AB|');
        fprintf('\na) for A = (1.2, 3.5), B = (12, 15): unitvec = ');disp(thirteenth([1.2, 3.5],[12, 15]));
        fprintf('\nb) for A = (-10, -4, 2.5), B = (-13, 6, -5): unitvec = ');disp(thirteenth([-10, -4, 2.5],[-13, 6, -5]));

        fprintf('\n14. crosspro=ABxCD');
        fprintf('\na) for a = (3, 11), b = (14, -7.3): crosspro = ');disp(fourteenth([3 11],[14 -7.3]));
        fprintf('\nb) for a = (-6, 14.2, 3), b = (6.3, -8, -5.6): crosspro = ');disp(fourteenth([-6 14.2 3],[6.3 -8 -5.6]));

        fprintf('\n15. S=0.5|ABxAC|');
        fprintf('\na) for A = (1, 2), B = (10, 3), C = (6, 11): S = %.2f', fiveteenth([1, 2], [10, 3], [6, 11]));
        fprintf('\nb) for A = (-1.5, -4.2, -3), B = (-5.1, 5.3, 2), C = (12.1, 0, -0.5): S = %.2f\n', fiveteenth([-1.5, -4.2, -3], [-5.1, 5.3, 2], [12.1, 0, -0.5]));

        fprintf('\n16. R=abc/4S;  L=2piR;');
        [L1,~]=sixteenth([1, 2], [10, 3], [6, 11]);
        fprintf('\na) for A = (1, 2), B = (10, 3), C = (6, 11): L = %.2f', L1);
        [L2,~]=sixteenth([-1.5, -4.2, -3], [-5.1, 5.3, 2], [12.1, 0, -0.5]);
        fprintf('\nb) for A = (-1.5, -4.2, -3), B = (-5.1, 5.3, 2), L = (12.1, 0, -0.5): S = %.2f\n', L2);
        clear L1 L2 R1 R2

        figure('Units', 'Normalized', 'OuterPosition', [1/3 0 1/3 1]);
        subplot(2, 1, 1)
        seventeenth([7.2, -2.9], [-1.8, 0.5]);
        title('17.1','Interpreter', 'latex');
        subplot(2, 1, 2)
        seventeenth([-0.9, -3.3], [0, 10]);
        title('17.2','Interpreter', 'latex');

        format longE
        fprintf('\n17. dex2bin\n');
        eighteenth(100);
        eighteenth(1002);
        eighteenth(52601);
        eighteenth(200090);

        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        ninetennth([1.5 3] , [9 10.5], [6 -3.8]);
        title('19','Interpreter', 'latex');

        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        subplot(2, 1, 1)
        twentyth(3.5, 2, 8.5, 3);
        title('20.1','Interpreter', 'latex');
        subplot(2, 1, 2)
        twentyth(-5, 1.5, 4, 8);
        title('20.2','Interpreter', 'latex');

        fprintf('\n21. r=r1+r2');
        [r,th]=twfirst(5, 23, 12, 40);
        fprintf('\na) for r1=5, th1=23° and r2=12, th2=40°: r=%0.1f th=%.1f°', [r,th]);
        [r,th]=twfirst(6, 80, 15, 125);
        fprintf('\nb) for r1=6, th1=80° and r2=15, th2=125°: r=%.0f th=%.0f°\n', [r,th]);
        clear r th

        fprintf('\n22. All the primes from %g to %g\n', [12, 80]);
        twsecond(12, 80)
        fprintf('\n    All the primes from %g to %g\n', [21, 63.5]);
        twsecond(21, 63.5)
        fprintf('\n    All the primes from %g to %g\n', [100, 200]);
        twsecond(100, 200)
        fprintf('\n    All the primes from %g to %g\n', [90, 50]);
        twsecond(90, 50)

        fprintf('\n23. GM for given array is %g\n', twthird([1.076, 1.113, 1.135, 1.103, 1.062, 1.032, 1.043, 1.036, 1.019, 1.036]));

        fprintf('\n24. CArt2Polar');
        x=14; y=9;
        [r,th]=twfourth(x, y);
        fprintf('\n for x=%.1f and y=%.1f: r=%0.1f th=%.1f°', [x,y,th,r]);
        x=-11; y=-20;
        [r,th]=twfourth(x, y);
        fprintf('\n for x=%.1f and y=%.1f: r=%0.1f th=%.1f°', [x,y,th,r]);
        x=-15; y=4;
        [r,th]=twfourth(x, y);
        fprintf('\n for x=%.1f and y=%.1f: r=%0.1f th=%.1f°', [x,y,th,r]);
        x=13.5; y=-23.5;
        [r,th]=twfourth(x, y);
        fprintf('\n for x=%.1f and y=%.1f: r=%0.1f th=%.1f°', [x,y,th,r]);
        clear x y
    
    elseif (Chapter == 8)
        p=[0.1 -0.2 -1 5 -41.5 235];
        x=-6:0.01:6;
        y=polyval(p,x);
        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(x,y);
        ylabel('$y$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('1','Interpreter', 'latex');
        grid on
        clear x y p

        p=[0.008 0 -1.8 -5.4 54];
        x=linspace(-14, 16);
        y=polyval(p,x);
        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(x,y);
        ylabel('$y$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('2','Interpreter', 'latex');
        grid on
        clear x y p

        a=[-3 0 5 -1];
        b=[1 2 -16 5];
        p=conv(a,b);
        fprintf('\n3.a*b=p\na = ');
        disp(a);
        fprintf('\nb = ');
        disp(b);
        fprintf('\np = ');
        disp(p);
        clear a b p

        a=[0 1.7 -0.5 0.7 -1.5];
        p=poly(a);
        x=linspace(-1.6, 1.8);
        y=polyval(p,x);
        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(x,y);
        ylabel('$y$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('4','Interpreter', 'latex');
        grid on
        clear a x y p

        a=[-10 -20 9 10 8 11 -3];
        b=[2 4 1];
        [p,q]=deconv(a,b);
        fprintf('\n5.a/b=p(q mod)\na = ');
        disp(a);
        fprintf('\nb = ');
        disp(b);
        fprintf('\np = ');
        disp(p);
        fprintf('\nq = ');
        disp(q);
        clear a b p q

        a=[-0.24 1.6 1.5 -7.41 -1.8 -4 -75.2 -91];
        b=[-0.8 0 5 6.5];
        [p,q]=deconv(a,b);
        fprintf('\n6.a/b=p(q mod)\na = ');
        disp(a);
        fprintf('\nb = ');
        disp(b);
        fprintf('\np = ');
        disp(p);
        fprintf('\nq = ');
        disp(q);
        clear a b p q

        r=[0 1];
        p=poly(r);
        res=6972;
        polyeq=[0 0 -res]+p;
        x=roots(polyeq);
        fprintf('\n7.x(x-1)=6972\n');
        fprintf('    x_%1i=%2.0f (x-1)_%1i=%2.0f\n',[(r+1); x'; (r+1); (x-1)']);
        clear r p res x polyeq

        r=[0 5 10];
        p=poly(r);
        res=10098;
        polyeq=[0 0 0 -res]+p;
        x=roots(polyeq);
        for i=1:length(x)
           if(abs(x(i))~=real(x(i)))
             x(i)=0;
           end
        end
        n=length(x);
        for i=1:n
           if(x(n-i+1)==0)
             x(n-i+1)=[];
           end
        end
        n=1:length(x);
        fprintf('\n8.x(x-5)(x-10)=10098\n');
        fprintf('    x_%1i=%2.0f (x-5)_%1i=%2.0f (x-10)_%1i=%2.0f\n',[n; x'; n; (x-5)'; n; (x-10)']);
        clear r p res x polyeq

        W=12212; m=0.284;
        V=W/m;
        x=[-1 240]; y=[-1 120]; z=[-2 80];
        Vol=conv(x,conv(y,z));
        polyeq=[0 0 0 (V-240*120*80)]+Vol;
        t=roots(polyeq);
         for i=1:length(t)
           if(abs(t(i))~=real(t(i)))
             t(i)=0;
           end
        end
        n=length(t);
        for i=1:n
           if(t(n-i+1)==0)
             t(n-i+1)=[];
           end
        end
        n=1:length(t);
        fprintf('\n9.Val=W/m Val=Vall-Vempt\n');
        fprintf('   t_%1i=%2.2f\n',[n; t']);
        clear i m n polyeq t V Vol W x y z

        W=42.27; m=0.101; R=10; h=24;
        V=W/m;
        a=[-1.5 R]; b=[-1 R];
        polyeq=4/3*conv(a,conv(a,a))+[0 h*conv(b,b)]-[0 0 0 V/pi];
        t=roots(polyeq);
        for i=1:length(t)
           if(abs(t(i))~=real(t(i)))
             t(i)=0;
           end
        end
        n=length(t);
        for i=1:n
           if(t(n-i+1)==0)
             t(n-i+1)=[];
           end
        end
        n=1:length(t);
        fprintf('\n10.Val=W/m Val=Vall-Vempt\n');
        fprintf('   t_%1i=%f\n',[n; t']);
        clear i m n polyeq t V Vol W a b h R

        L=20*12;
        a=[1 0]; b=[1 15]; c=([0 L]-4*a-4*b)/4;
        V=conv(a,conv(b,c));
        x0=roots(polyder(V));
        n=length(x0);
        for i=1:length(x0)
           if(x0(n-i+1)<0)
                x0(n-i+1)=[];
           end
        end
        V0=polyval(V,x0);
        x=linspace(0,45/2);    
        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        y=polyval(V,x);
        plot(x,y);
        ylabel('$Volume=V(x)$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('11','Interpreter', 'latex');
        hold on
        plot(x0, V0, 'sb', 'MarkerSize', 7,'MarkerFaceColor', 'r')
        text(x0-5,V0+200,['Max Value = ', num2str(V0)]);
        ylim([y(1),V0+500]);
        grid on
        clear L a b c V x0 x V0 y i

        a=[-2 40]; b=[-2 22]; c=[1 0];
        V=conv(a,conv(b,c));
        x=linspace(0, 11);
        x1=roots(polyder(V));
        n=length(x1);
        for i=1:length(x1)
           if(polyval(V,x1(n-i+1))<0)
                x1(n-i+1)=[];
           end
        end
        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        y=polyval(V,x);
        plot(x,y);
        ylabel('$Volume=V(x)$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('12','Interpreter', 'latex');
        hold on
        x0=roots(V-[0 0 0 1000]);
        n=length(x0);
        for i=1:length(x0)
           if(x0(n-i+1)<0||x0(n-i+1)>11)
                x0(n-i+1)=[];
           end
        end
        plot(x0, polyval(V,x0), 'sb', 'MarkerSize', 7,'MarkerFaceColor', 'r')
        text(x0*0.8,polyval(V,x0)*0.9, 'V=1000');
        plot(x1, polyval(V,x1), 'ob', 'MarkerSize', 7,'MarkerFaceColor', 'r')
        text(x1*0.8,polyval(V,x1)/0.95, ['Max Value = ', num2str(polyval(V,x1))]);
        grid on
        clear a b c i n V x x0 x1 y

        F1=[2 0 -3 -9 11 -8 4];
        F2=[5 0 7 -10];
        FA=polyadd(F1, F2, 'add');
        FS=polyadd(F1, F2, 'sub');
        fprintf('\n13.Sum and diff of polynoms\n');
        fprintf('   F1+F2 = ');
        disp(FA);
        fprintf('\n   F1-F2 = ');
        disp(FS);
        clear FA FS

        FMLT=conv(F1,F2);
        FM=polymult(F1, F2);
        fprintf('\n14.Mult of polynoms\n');
        fprintf('   F1*F2   =   ');
        disp(FM);
        fprintf('\n   conv(F1,F2)=');
        disp(FMLT);
        clear F1 F2 FM FMLT

        F1=[3 -7 14]; F2=[-5 -11 15];
        [x,y,w]=maxormin(F1);
        fprintf('\n15.x0=%6.4g y0=%6.4g w=%2g\n', [x,y,w]);
        [x,y,w]=maxormin(F2);
        fprintf('   x0=%6.4g y0=%6.4g w=%2g\n', [x,y,w]);
        clear F1 F2 w x y

        R=9; 
        V=pi/3*conv([-1 0 R^2],[1 R]);
        h=linspace(-R,R);
        Vol=polyval(V,h);
        H=roots(V-[0 0 0 500]);
        h0=roots(polyder(V));
        n=length(H);
        for i=1:length(H)
           if(H(n-i+1)<-R||H(n-i+1)>R)
                H(n-i+1)=[];
           end
        end
        n=length(H);
        for i=1:length(h0)
           if(polyval(polyder(polyder(V)),h0(n-i+1))>0)
                h0(n-i+1)=[];
           end
        end
        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(h,Vol);
        xlim([-R, R]);
        ylabel('$Volume=V(h)$','Interpreter', 'latex');
        xlabel('$h$','Interpreter', 'latex');
        title('16','Interpreter', 'latex');
        grid on
        hold on
        plot(H, polyval(V,H), 'sb', 'MarkerSize', 7,'MarkerFaceColor', 'r')
        text(H(1)-2, polyval(V,H(1))-50, ['Value = ', num2str(polyval(V,H(1)))]);
        text(H(2)-2, polyval(V,H(2))-50, ['Value = ', num2str(polyval(V,H(2)))]);
        plot(h0, polyval(V,h0), 'ob', 'MarkerSize', 7,'MarkerFaceColor', 'r')
        text(h0-2, polyval(V,h0)-50, ['Value = ', num2str(polyval(V,h0))]);
        clear h H h0 i n R V Vol

        P=[3 5.5]; d0=28;
        Y=conv([1 -3],[1 -3])*1.5+[0 0 1];
        x=linspace(-2,8);
        y=polyval(Y,x);
        D=[0 0 conv([-1 P(1)],[-1 P(1)])]+...
            conv((-Y+[0 0 P(2)]),(-Y+[0 0 P(2)]));
        d=sqrt(polyval(D,x));
        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(x,d, 'LineWidth', 2);
        ylabel('$d=dist(P,Q)$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('17','Interpreter', 'latex');
        grid on
        hold on
        plot(x, y, '--b');
        x0=roots(D-[0 0 0 0 d0^2]);
        n=length(x0);
        for i=1:n
           if(~isreal(x0(n-i+1))==1)
             x0(n-i+1)=[];
           end
        end
        text(x0(1)-1, sqrt(polyval(D,(x0(1))))-2, ['Value = ', num2str(sqrt(polyval(D,(x0(1)))))]);
        text(x0(2)-1, sqrt(polyval(D,(x0(2))))-2, ['Value = ', num2str(sqrt(polyval(D,(x0(2)))))]);
        plot(x0, sqrt(polyval(D,(x0))), 'ob', 'MarkerSize', 7,'MarkerFaceColor', 'r')
        plot(x,28*x.^0);
        xm=roots(polyder(D));
        n=length(xm);
        for i=1:n
           if(polyval(polyder(polyder(D)),xm(n-i+1))<0)
                xm(n-i+1)=[];
           end
        end
        plot(xm, sqrt(polyval(D,(xm))), 'sb', 'MarkerSize', 7,'MarkerFaceColor', 'r')
        text(xm(1)-1, sqrt(polyval(D,(xm(1))))-1, ['minValue = ', num2str(sqrt(polyval(D,(xm(1)))))]);
        text(xm(2)-1, sqrt(polyval(D,(xm(2))))-1, ['minValue = ', num2str(sqrt(polyval(D,(xm(2)))))]);
        xlim([x(1)-1, x(100)+1])
        clear d D d0 i n P x x0 xm y Y

        x=[2 5  6  8  9 13 15];
        y=[7 8 10 11 12 14 15];
        p=polyfit(x,y,1);
        X=linspace(x(1),x(length(x)));
        Y=polyval(p,X);
        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(X,Y, x,y, 'sk');
        ylabel('$y$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('18','Interpreter', 'latex');
        grid on
        clear x y X Y p

        h0=5000;
        h=[0 600 1500 2300 3000 6100 7900];
        t=[100 98.8 95.1 92.2 90 81.2 75.6];
        p=polyfit(h,t,1);
        H=linspace(h(1),h(length(h)));
        T=polyval(p,H);
        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(H,T, h,t, 'sk', h0,polyval(p, h0), 'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        ylabel('$T$','Interpreter', 'latex');
        xlabel('$h$','Interpreter', 'latex');
        text(h0+100,polyval(p, h0)+1,['T(5000)=', num2str(polyval(p, h0))]);
        title('19','Interpreter', 'latex');
        grid on
        clear h t H T p h0

        t0=1915;
        t=[1815 1845 1875 1905 1935 1965];
        p=[8.3 19.7 44.4 83.2 127.1 190.9];
        n=polyfit(t,p,2);
        T=linspace(t(1),t(length(t)));
        P=polyval(n,T);
        figure('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(T,P, t,p, 'sk', t0,polyval(n, t0), 'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        ylabel('$Population/10^6$','Interpreter', 'latex');
        xlabel('$Years$','Interpreter', 'latex');
        text(t0+15,polyval(n, t0),['T(5000)=', num2str(polyval(n, t0))]);
        title('20','Interpreter', 'latex');
        grid on
        clear t p n T P t0
    
    elseif(Chapter == 9)
        figure('Units', 'normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        x=linspace(0,15);
        y=exp(0.3*x)-x.^2;
        plot(x,y);
        ylabel('$y$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('$1. exp(0.3*x)-x^2=-4$','Interpreter', 'latex');
        grid on;
        x1=fzero('exp(0.3*x)-x^2+4',2);
        hold on;
        plot(x,x.^0*-4,'--');
        plot(x1,-4,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        text(x1,-4*1.1,['   (',num2str(x1),',',num2str(-4), ')']);
        clear x y x1

        figure('Units', 'normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        x=linspace(0,5);
        y=2*cos(x)-0.5*x;
        plot(x,y);
        ylabel('$y$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('$2. 2cos(x)-0.5*x=1$','Interpreter', 'latex');
        grid on;
        x1=fzero('2*cos(x)-0.5*x-1',2);
        hold on;
        plot(x,x.^0*1,'--');
        plot(x1,1,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        text(x1,1*1.1,['   (',num2str(x1),',',num2str(1), ')']);
        clear x x1 val


        figure('Units', 'normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        x=linspace(0,7);
        y=x.^3-5*x.^2.5+exp(0.9*x)+4*(x+1);
        val=-2;
        plot(x,y);
        ylabel('$y$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('$3. x^3-5x^{2.5}+exp(0.9x)+4(x+1)=-2$','Interpreter', 'latex');
        grid on;
        hold on;
        plot(x,x.^0*val,'--');
        x1=fzero('x.^3-5*x.^2.5+exp(0.9*x)+4*(x+1)+2',2);
        plot(x1,val,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        text(x1,val*1.1,['   (',num2str(x1),',',num2str(val), ')']);
        x2=fzero('x.^3-5*x.^2.5+exp(0.9*x)+4*(x+1)+2',7);
        plot(x2,val,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        text(x2,val*1.1,['   (',num2str(x2),',',num2str(val), ')']);
        ylim([-45,10]);
        clear x y x1 x2 val

        figure('Units', 'normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        x=linspace(1,4);
        y=x.^2-5*x.*sin(3*x)+3;
        val=0;
        plot(x,y);
        ylabel('$y$','Interpreter', 'latex');
        xlabel('$x$','Interpreter', 'latex');
        title('$4. x^2-5xsin(3x)+3=0$','Interpreter', 'latex');
        grid on;
        hold on;
        plot(x,x.^0*val,'--');
        x1=fzero('x.^2-5*x.*sin(3*x)+3',2);
        plot(x1,val,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        text(x1,val*1.1,['   (',num2str(x1),',',num2str(val), ')']);
        x2=fzero('x.^2-5*x.*sin(3*x)+3',3);
        plot(x2,val,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        text(x2,val*1.1,['   (',num2str(x2),',',num2str(val), ')']);
        clear x y x1 x2 val

        m=25; g=9.81; mu=0.55; F=150;
        Tet=fzero(@(tet)(mu*m*g)/(cosd(tet)+mu*sind(tet))-F,45);
        fprintf('\n5.For F=%g, tet=%3g°\n', [F Tet]);
        clear F g m mu Tet

        a=0.22; b=0.08; K1=1600; K2=1e5;
        x=linspace(0,0.25);
        L=@(x)sqrt(a^2+(b+x).^2);
        L0=sqrt(a^2+b^2);
        u=@(x) L(x)-L0;
        Fs=@(x)K1*u(x)+K2*u(x).^3;
        W=2*Fs(x).*(b+x)./L(x);
        figure('Units', 'normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(x,W);
        ylabel('$W,Newtons$','Interpreter', 'latex');
        xlabel('$x, meters$','Interpreter', 'latex');
        title('$6. W=2F_s{(b+x)\over L}$','Interpreter', 'latex');
        grid on;
        hold on;
        x0=fzero(@(x)2*Fs(x).*(b+x)./L(x)-400,0.175);
        plot(x0,400,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        text(x0,400*1.05,['   (',num2str(x0),',',num2str(400), ')']);
        clear a b K1 K2 x L L0 u Fs W x0

        g=9.81; V=0.8; M=0.1; ro=1000; beta=10; teta=10; C=1;
        f=@(d)sqrt(16*M*g/(pi*C*ro*d^2))/sqrt(1-8*M*tand(beta)^2/(pi*d^3*C*ro*sind(teta)))-V;
        diam=fzero(f, 0.1);
        fprintf('\n7.For V=%g, d=%3.3g\n', [V diam]);
        clear g V M ro beta teta C f diam

        Is=1e-12; k=1.38e-23; q=1.6e-19; T=297; vs=2; R=1000;
        I=@(wd)Is*(exp(q*wd/(k*T))-1);
        II=@(wd)(vs-wd)/R;
    %     W=linspace(0.5,0.6);
    %     I=I(W);
    %     II=II(W);
    %     figure('Units', 'normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
    %     plot(W,I,W,II);
        w=fzero(@(wd)I(wd)-II(wd),1);
        fprintf('\n8.For v_s=%g, w_d=%3.3g\n', [vs w]);
        clear Is k q T vs R I II w

        f=@(x,s) s*3*(x-1/4)./(1+7/2*(4/5*x-3/10).^2);
        x=linspace(-5,5);
        y=f(x,1);
        figure('Units', 'normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(x,y);
        ylabel('$W,Newtons$','Interpreter', 'latex');
        xlabel('$x, meters$','Interpreter', 'latex');
        title('$10. MinAndMaxOf{3(x-0.25)\over 1+3.5(0.8x-0.3)^2}$','Interpreter', 'latex');
        grid on;
        hold on;
        [x1,y1]=fminbnd(@(x)f(x,1),-2,2);
        plot(x1,y1,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        text(x1,y1,['   min value(',num2str(x1),',',num2str(y1), ')']);
        [x2,y2]=fminbnd(@(x)f(x,-1),-2,2);
        plot(x2,-y2,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'r');
        text(x2,-y2,['   max value(',num2str(x2),',',num2str(-y2), ')']);
        clear f x y x1 y1 x2 y2

        R2=@(R1) 2*R1;
        V=250;
        h=@(R1) 3*V/(pi*(R1^2+R2(R1)^2+R1*R2(R1)));
        S=@(R1) pi*(R1+R2(R1))*sqrt((R1-R2(R1))^2+h(R1)^2)+pi*(R1^2+R2(R1)^2);
        [R1,Smin]=fminbnd(S, 0, 100);
        hmin=h(R1);
        R2min=R2(R1);
        fprintf('\n10.For V=%g: S_min=%3.2g, h=%3.2g R1=%3.2g R2=%3.2g\n', ...
            [V Smin hmin R1 R2min]);
        clear R2 V h S R1 Smin hmin R2min

        m=25; g=9.81; mu=0.55; 
        F=@(tet)(mu*m*g)/(cosd(tet)+mu*sind(tet));
        [tet, Fmin]=fminbnd(@(tet) F(tet), 0, 90);
        fprintf('\n11.For Fmin=%g, tet=%3g°\n', [Fmin tet]);
        clear m g mu F tet Fmin

        R=14;
        h=@(r)sqrt(R^2-r^2);
        V=@(r) -pi*r^2*2*h(r);
        [r,V0]=fminbnd(V, 0 , R);
        fprintf('\n12.For Sphere with R=%g, Vcilin=%2g h=%2g r=%2g\n', [R -V0 2*h(r) r]);
        clear R h V r V0

        y=@(x) 5*sqrt(1-x^2/19^2);
        S=@(x) 2*y(x)*2*x;
        [x0,S0]=fminbnd(@(x) -S(x),0,19);
        y0=y(x0);
        a=2*x0;
        b=2*y0;
        fprintf('\n13.For ellipse x^2/19^2+y^2/5^2=1 Smax=%g, a=%g, b=%g\n', [-S0 a b]);
        clear a b x0 y0 x y S S0

        c=3e8; h=6.63e-34; k=1.38e-23; T=1500;
        lmbd=linspace(0.2e-6,6e-6);
        R=@(l)2*pi*c^2*h./l.^5./(exp(h*c./(l*k*T))-1);
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(lmbd,R(lmbd));
        ylabel('$R$','Interpreter', 'latex');
        xlabel('$\lambda, meters$','Interpreter', 'latex');
        title('$14. R={2pic^2h\over \lambda^5}{1\over e^{(h*c/(\lambda*k*T))}-1}$','Interpreter', 'latex');
        grid on;
        hold on;
        n=10^6;
        [l1, ~]=fminbnd(@(l)-R(l*n),0.2e-6*n,6e-6*n);
        plot(l1/n*10,R(l1/n*10),'ok', 'MarkerSize', 7,'MarkerFaceColor', 'b');
        text(l1/n*10,R(l1/n*10),'   Maxvalue');
        clear c h k l1 lmbd n R R1 T

        L=108; Lc=68; W=250;
        T=@(d) W*L*Lc./(d.*sqrt(Lc^2-d.^2));
        d=linspace(0,Lc);
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(d,T(d));
        ylabel('$T, newtons$','Interpreter', 'latex');
        xlabel('$d, inches$','Interpreter', 'latex');
        title('$15. T={WLL_c\over d\sqrt{Lc^2-d^2}}$','Interpreter', 'latex');
        grid on;
        hold on;
        [d0,T0]=fminbnd(T,0,Lc);
        plot(d0,T0,'ok', 'MarkerSize', 7,'MarkerFaceColor', 'b');
        text(d0,T0*2.5,'   minvalue');
        xlim([d(1),d(100)]);
        ylim([0,5e3]);
        clear d d0 L Lc T T0 W

        f=@(x)0.5*x./(1+2*sqrt(x));
        I=integral(f,2,10);
        fprintf('\n16. a) Integral 0.5x/(1+2sqrt(x)) from 2 to 10 =%3g', I);
        f=@(x)0.5+cos(1.2*x)./(x+2).^2;
        I=integral(f,0,9);
        fprintf('\n    b) Integral 0.5+cos(1.2x)/(x+2)^2 from 0 to 9=%3g\n', I);
        clear f I

        f=@(x)exp(x)./x.^3;
        I=integral(f,1,8);
        fprintf('\n17. a) Integral e^x/x^3 from 1 to 8 =%3g', I);
        f=@(x)cos(x).*exp(sqrt(x));
        I=integral(f,0,4*pi);
        fprintf('\n    b) Integral cos(x)*exp(sqrt(x)) from 0 to 4pi =%3g\n', I);
        clear f I

        t=0:7;
        v=[0 14 39 69 95 114 129 139];
        I=trapz(t(1:7),v(1:7));
        fprintf('\n18. Distante beaten by car=%3g\n', I);
        clear t v I

        syms f(x)
        f(x)=@(x)693.9-68.8*cosh(x/99.7);
        F=diff(f,x);
        x=linspace(-299.26, 299.25);
        F=F(x);
        L=trapz(x, sqrt(1+F.^2));
        fprintf('\n19. Length of arc = %3g\n', L);
        clear f F x L

        R=0.25; n=7; vmax=80;
        v=@(r)vmax*(1-r/R).^(1/n);
        Q=integral(@(r) 2*pi*v(r).*r,0,R);
        fprintf('\n20. consumption of water = %3g\n', Q);
        clear R n vmax v Q


        sg=300e-6; eps=8.85e-12; z=0.05; R=0.06;
        E=sg*z/4/eps*integral(@(r)((z^2+r.^2).^(-3/2)).*r*2,0,R);
        fprintf('\n21. E(z) = %f\n', E);
        clear sg eps z R E

        t=linspace(0,2*pi);
        b=5;
        x=@(t)2*b*cos(t)-b*cos(2*t);
        y=@(t)2*b*sin(t)-b*sin(2*t);
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot(x(t),y(t));
        title('$22. Cardioida$','Interpreter', 'latex');
        grid on;
        syms x(t)
        T=linspace(0,2*pi);
        x(t)=@(t)2*b*cos(t)-b*cos(2*t);
        X=diff(x,t);
        X=X(T); Y=y(T);
        L=trapz(T,sqrt(X.^2+Y.^2));
        text(-2*b,0,['Length equals to ', num2str(round(double(L)*100)/100)],'Interpreter', 'latex');
        xlim([-20,10]);
        ylim([-15,15]);
        clear t b x y T X Y L

        g0=9.81; R=6371e3;
        g=@(y)R^2./(R+y).^2*g0;
        m=500; h=800e3;
        U=m*integral(g,0,h);
        fprintf('\n23. Required energy = %3g\n', U);
        clear g0 R g m h U

        h=[0 40 96 140 147 121 117 139 140 62 18 0];
        dx=40;
        n=length(h);
        x=0:40:40*(n-1);
        S=trapz(x,h);
        hmean=S/(n-1)/40;
        fprintf('\n24. Approximate square = %g, h_mean = %g', [S hmean]);
    
    elseif(Chapter == 10)
        t=linspace(0,30, 1000);
        x=((t-15)/100+1).*sin(3*t);
        y=((t-15)/100+1).*cos(0.8*t);
        z=0.4*t.^1.5;
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot3(x,y,z);
        title('$1$','Interpreter', 'latex');
        xlabel('x','Interpreter', 'latex');
        ylabel('y','Interpreter', 'latex');
        zlabel('z','Interpreter', 'latex');
        grid on;
        axis equal
        xlim([-35, 35]);
        ylim([-35, 35]);
        zlim([0, 70]);
        clear x y z t

        a=20; b=10; h=18; n=3;
        t=linspace(0,n*2*pi, 1000);
        r=a*b./sqrt((b*cos(t)).^2+(a*sin(t)).^2);
        x=r.*cos(t);
        y=r.*sin(t);
        z=h*t/pi/n;
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot3(x,y,z);
        title('$2$','Interpreter', 'latex');
        xlabel('x','Interpreter', 'latex');
        ylabel('y','Interpreter', 'latex');
        zlabel('z','Interpreter', 'latex');
        grid on;
        axis equal
        xlim([-20, 20]);
        ylim([-20, 20]);
        zlim([0, 40]);
        clear x y z t n r a b h

        vphi=5*pi/180; vtetha=8*pi/180; vr=0.6;
        t=linspace(0,10,500);
        phi=t*vphi;
        tetha=t*vtetha;
        r=t*vr;
        [x,y,z]=sph2cart(tetha, phi,r);
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        plot3(x,y,z);
        title('$3$','Interpreter', 'latex');
        xlabel('x','Interpreter', 'latex');
        ylabel('y','Interpreter', 'latex');
        zlabel('z','Interpreter', 'latex');
        grid on;
        axis equal
        xlim([-4, 4]);
        ylim([-4, 4]);
        zlim([0, 8]);
        clear x y z phi tetha r vphi vtetha vr t

        x=linspace(-3,3);
        y=linspace(-3,3);
        [X,Y]=meshgrid(x,y);
        Z=Y.^2/4-2*sin(1.5*X);
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        surf(X,Y,Z, 'EdgeColor','none');
        title('$4$','Interpreter', 'latex');
        xlabel('x','Interpreter', 'latex');
        ylabel('y','Interpreter', 'latex');
        zlabel('z','Interpreter', 'latex');
        grid on;
        axis equal
        xlim([x(1), x(100)]);
        ylim([y(1), y(100)]);
        zlim([-3, 5]);
        colormap('jet');
        clear x y z X Y Z

        x=linspace(-2,2);
        y=linspace(-2,2);
        [X,Y]=meshgrid(x,y);
        Z=0.5*X.^2+0.5*Y.^2;
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        surf(X,Y,Z, 'EdgeColor','none');
        title('$5$','Interpreter', 'latex');
        xlabel('x','Interpreter', 'latex');
        ylabel('y','Interpreter', 'latex');
        zlabel('z','Interpreter', 'latex');
        grid on;
        axis equal
        xlim([x(1), x(100)]);
        ylim([y(1), y(100)]);
        zlim([0, 4]);
        colormap('bone');
        clear x y Z X Y

        x=linspace(-5,5);
        y=linspace(-5,5);
        [X,Y]=meshgrid(x,y);
        Z=-cos(2*sqrt(X.^2+Y.^2))./exp(0.2*sqrt(X.^2+Y.^2));
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        surf(X,Y,Z);
        title('$6$','Interpreter', 'latex');
        xlabel('x','Interpreter', 'latex');
        ylabel('y','Interpreter', 'latex');
        zlabel('z','Interpreter', 'latex');
        grid on;
        axis equal
        xlim([x(1), x(100)]);
        ylim([y(1), y(100)]);
        zlim([-5, 5]);
        colormap('cool');
        clear x y Z X Y

        x=linspace(-2*pi,2*pi);
        y=linspace(-pi,pi);
        [X,Y]=meshgrid(x,y);
        Z=cos(X).*cos(sqrt(X.^2+Y.^2)).*exp(-abs(0.2*X));
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        surf(X,Y,Z);
        title('$7$','Interpreter', 'latex');
        xlabel('x','Interpreter', 'latex');
        ylabel('y','Interpreter', 'latex');
        zlabel('z','Interpreter', 'latex');
        grid on;
        axis equal
        xlim([x(1), x(100)]);
        ylim([y(1), y(100)]);
        zlim([-5, 5]);
        colormap('hsv');
        clear x y Z X Y

        [r,tet]=meshgrid(linspace(0,2),linspace(0,2*pi));
        x1=r.*cos(tet);
        y1=r.*sin(tet);
        z1=4*r;
        [phi,tet]=meshgrid(linspace(0,pi/2),linspace(0,2*pi));
        x2=2*cos(tet).*sin(phi);
        y2=2*sin(tet).*sin(phi);
        z2=8+2*cos(phi);
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        grid on;
        surf(x1,y1,z1,'EdgeColor', '#EDB120');
        hold on
        surf(x2,y2,z2,'EdgeColor', 'k');
        axis equal
        colormap('white');
        title('$8$','Interpreter', 'latex');
        xlabel('x','Interpreter', 'latex');
        ylabel('y','Interpreter', 'latex');
        zlabel('z','Interpreter', 'latex');
        clear r tet x1 y1 z1 phi tet x2 y2 z2 s1 s2

        R=0.08206; n=1.5; a=1.39; b=0.03913;
        [V,T]=meshgrid(linspace(0.3,1.2),linspace(273,473));
        P=n*R*T./(V-b*n)-n^2*a./V^2;
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        grid on;
        surf(V,T,P, 'EdgeColor','none');
        title('$9. Van-derh-Walts law$','Interpreter', 'latex');
        xlabel('$V, L$','Interpreter', 'latex');
        ylabel('$T, K$','Interpreter', 'latex');
        zlabel('$P, atm$','Interpreter', 'latex');
        %axis equal
        %colormap('white');
        colorbar;
        clear R n a b V T P

        R=8.31; M=0.032;
        [v,T]=meshgrid(linspace(0,1000),linspace(70,320));
        P=4*pi*(M./(2*pi*R*T)).^1.5.*v^2.*exp(-M*v.^2./(2*R*T));
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        grid on;
        surf(v,T,P,'EdgeColor','none');
        title('$10. Distribution of speed$','Interpreter', 'latex');
        xlabel('$v, m/s$','Interpreter', 'latex');
        ylabel('$T, K$','Interpreter', 'latex');
        zlabel('$P$','Interpreter', 'latex');
        colormap('parula');
        clear R M v T P

        C1=3.742e8; C2=1.439e4;
        [l,T]=meshgrid(linspace(0.1,10),linspace(10,2000));
        E=C1./(l.^5.*(exp(C2./(l.*T))-1));
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        grid on;
        surf(l,T,E,'EdgeColor','none');
        set(gca,'xscale','log');
        title('$11. PlanksLaw$','Interpreter', 'latex');
        xlabel('$\lambda, mkm*T$','Interpreter', 'latex');
        ylabel('$T, K$','Interpreter', 'latex');
        zlabel('$E, {W\over m^2mkm}$','Interpreter', 'latex');
        colormap('cool');
        clear C1 C2 l T E

        [w,d]=meshgrid(linspace(0,8),linspace(0,4));
        k=1; n=0.05; S=0.001;
        %Q=((w.*d./(w+2*d))^2).^(1/3);
        Q=k*d.*w/n.*(((w.*d./(w+2*d))^2).^(1/3))*sqrt(S);
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        grid on;
        surf(w,d,Q,'EdgeColor','none');
        title('$12. MannigsLaw$','Interpreter', 'latex');
        xlabel('$width of channel, m$','Interpreter', 'latex');
        ylabel('$depth of channel, m$','Interpreter', 'latex');
        zlabel('$Q$','Interpreter', 'latex');
        colormap('cool');
        clear w d k n S Q

        C=15e-6; L=240e-3; vm=24;
        [f,R]=meshgrid(linspace(60,110,500),linspace(10,40,500));
        w=2*pi*f;
    %   vs=vm*sin(w.*t);
        I=vm./sqrt(R.^2+(w.*L-1./(w*C)).^2);
        figure ('Units', 'Normalized', 'OuterPosition', [1/4 0 1/2 1]);
        grid on;
        subplot(2,1,1)
        surf(f,R,I,'EdgeColor','none');
        title('$13. a)I=I(f_d,R)$','Interpreter', 'latex');
        xlabel('$f_d, s^{-1}$','Interpreter', 'latex');
        ylabel('$R, Ohms$','Interpreter', 'latex');
        zlabel('$I, A$','Interpreter', 'latex');
        colormap('winter');

        subplot(2,1,2)
        plot(f(1,:),I(1,:,1))
        title('$13. b)I=I(f_d,10)$','Interpreter', 'latex');
        xlabel('$f_d, s^{-1}$','Interpreter', 'latex');
        ylabel('$I, A$','Interpreter', 'latex');

        v0=1/(2*pi*sqrt(L*C));
        w0=2*pi*v0; 
        R=10;
        I0=vm./sqrt(R^2+(w0*L-1/(w0*C))^2);
        hold on 
        plot(v0, I0, 'ok', 'MarkerSize', 7, 'MarkerFaceColor','r');
        text(65, 2, ['Calculated value = ', num2str(v0)],'Color','r','FontSize',14);

        [~,n]=max(I(1,:,1));
        plot(f(1,n),I(1,n,1), 'sk', 'MarkerSize', 7, 'MarkerFaceColor','b');
        text(93, 2, ['Graphs value = ', num2str(f(1,n))],'Color','b','FontSize',14);
    
    elseif(Chapter==11)%19
        syms x
        S1=x^2*(x-6)+4*(3*x-2);
        S2=(x+2)^2-8*x;
        V1=simplify(S1*S2);
        V2=simplify(S1/S2);
        V3=simplify(S1+S2);
        fprintf('\n1.S1 = ');
        disp(S1);
        fprintf('  S2 = ');
        disp(S2);
        fprintf('  S1*S2 = ');
        disp(V1);
        fprintf('  S1/S2 = ');
        disp(V2);
        fprintf('  S1+S2 = ');
        disp(V3);
        fprintf('  S1+S2(x=5) = ');
        disp(subs(V3,5));
        clear S1 S2 V1 V2 V3 V4 x

        syms x
        S1=x*(x^2+6*x+12)+8;
        S2=(x-3)^2+10*x-5;
        V1=simplify(S1*S2);
        V2=simplify(S1/S2);
        V3=simplify(S1+S2);
        fprintf('\n2.S1 = ');
        disp(S1);
        fprintf('  S2 = ');
        disp(S2);
        fprintf('  S1*S2 = ');
        disp(V1);
        fprintf('  S1/S2 = ');
        disp(V2);
        fprintf('  S1+S2 = ');
        disp(V3);
        fprintf('  S1+S2(x=3) = ');
        disp(subs(V3,3));
        clear S1 S2 V1 V2 V3 V4 x

        syms x y
        S=x+sqrt(x)*y^2+y^4;
        T=sqrt(x)-y^2;
        V=simplify(S*T);
        fprintf('\n3.S = ');
        disp(S);
        fprintf('  T = ');
        disp(T);
        fprintf('  S*T = ');
        disp(V);
        fprintf('  S*T(y=2) = ');
        disp(subs(V,y,2));
        clear x y S V V1

        syms x
        t=[-2 -0.5 2 4.5];
        S=1;
        for i=1:length(t)
            S=S*(x-t(i)); 
        end
        S=expand(S);
        fprintf('\n4.a)S = ');
        disp(S);

        f=x^6-6.5*x^5-58*x^4+167.5*x^3+728*x^2-890*x-1400;
        f=factor(f);
        L(1:length(f)-1)=0;
        for i =2:length(f)
            L(i-1)=root(f(i));
        end
        i=1:length(f)-1;
        fprintf('4.b)For x^6-6.5*x^5-58*x^4+167.5*x^3+728*x^2-890*x-1400\n');
        fprintf('x%i = %g;  ',[i' L']);
        %disp(expand(f));
        clear f i S t T x L

        syms x y
        f=sin(4*x);
        nof=4*sin(x)*cos(x)-8*sin(x)^3*cos(x);
        sol=simplify(f-nof);
        fprintf('\n\n5.a)sin(4*x) = ');
        disp(expand(f));
        fprintf('    4*sin(x)*cos(x)-8*sin(x)^3*cos(x) = ');
        disp(simplify(4*sin(x)*cos(x)-8*sin(x)^3*cos(x)));
        fprintf('    sin(4*x)-4*sin(x)*cos(x)-8*sin(x)^3*cos(x) = ');
        disp(sol);

        j=cos(x)*cos(y);
        noj=0.5*(cos(x-y)+cos(x+y));
        sol=simplify(j-noj);
        fprintf('\n5.b)cos(x)*cos(y) = ');
        disp(expand(j));
        fprintf('    0.5*(cos(x-y)+cos(x+y)) = ');
        disp(simplify(0.5*(cos(x-y)+cos(x+y))));
        fprintf('    cos(x)*cos(y)-0.5*(cos(x-y)+cos(x+y)) = ');
        disp(sol);
        clear x y f j nof noj

        syms x y z
        f=tan(3*x);
        nof=(3*tan(x)-tan(x)^3)/(1-3*tan(x)^2);
        sol=simplify(f-nof);
        fprintf('\n6.a)tan(3x) = ');
        disp(expand(f));
        fprintf('    (3*tan(x)-tan(x)^3)/(1-3*tan(x)^2) = ');
        disp(simplify((3*tan(x)-tan(x)^3)/(1-3*tan(x)^2)));
        fprintf('    tan(3x)-(3*tan(x)-tan(x)^3)/(1-3*tan(x)^2) = ');
        disp(sol);


        j=sin(x+y+z);
        noj=sin(x)*cos(y)*cos(z)+sin(z)*cos(x)*cos(y)+...
            sin(y)*cos(z)*cos(x)-sin(x)*sin(y)*sin(z);
        sol=simplify(j-noj);
        fprintf('\n6.b)sin(x+y+z) = ');
        disp(expand(j));
        fprintf('    sin(x)*cos(y)*cos(z)+sin(z)*cos(x)*cos(y)+sin(y)*cos(z)*cos(x)-sin(x)*sin(y)*sin(z) = ');
        disp(simplify(sin(x)*cos(y)*cos(z)+sin(z)*cos(x)*cos(y)+...
            sin(y)*cos(z)*cos(x)-sin(x)*sin(y)*sin(z)));
        fprintf('    sin(x+y+z) - sin(x)*cos(y)*cos(z)+sin(z)*cos(x)*cos(y)+sin(y)*cos(z)*cos(x)-sin(x)*sin(y)*sin(z) = ');
        disp(sol);
        clear x y z f nof j noj 

        syms x y t
        x=3*t/(1+t^3);
        y=3*t^2/(1+t^3);
        left=x^3+y^3;
        right=3*x*y;
        sol=simplify(left-right);
        fprintf('\n7.x^3+y^3=3xy\n  x^3 = ');
        disp(x^3);
        fprintf('  y^3 = ');
        disp(y^3);
        fprintf('  x^3+y^3 = ');
        disp(left);
        fprintf('  x^3+y^3 = ');
        disp(simplify(left));
        fprintf('  3xy = ');
        disp(right);
        fprintf('  x^3+y^3-3xy= ');
        disp(sol);
        clear x y t
        syms x y
        left=x^3+y^3;
        right=3*x*y;
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        fimplicit(left-right,[-2 2 -3 2])
        grid on
        xlabel('x', 'Interpreter', 'latex')
        ylabel('y', 'Interpreter', 'latex');
        title('7. $Curve$ $of$ $Descartes$', 'Interpreter', 'latex')
        clear left right x y sol

        %h=10; V=1050;
        fprintf('\n8.Volume=pi*R^2*h+2/3*pi*R^3\n   r=');
        syms R h V real
        Vol=pi*R^2*h+2/3*pi*R^3;
        [r,~,~]=solve(Vol==V,R, 'Real', true,'ReturnConditions',true);
        r=r(1);
        disp(r);
        r=double(subs(r,{h, V},{10,1050}));
        fprintf('   r =');
        disp(r);
        fprintf('   at r=%3.5g Volume=pi*R^2*h+2/3*pi*R^3 = %3.3f\n', [r,double(subs(Vol,{R,h},{r,10}))]);

        clear h V r Vol R

        syms T a v b T0 vm
        F=(T+a)*(v+b)-(T0+a)*b;
        fprintf('\n9.(T+a)*(v+b)=(T0+a)b\n  Vmax = ');
        Vm=solve(subs(F,T,0),v);
        disp(Vm);
        B=solve(Vm==vm,b);
        F=simplify(solve(subs(F,b,B),v));
        fprintf('  T in terms w/o b: T = ');
        disp(F);
        clear T a v b T0 vm F B Vm

        syms x y
        E1=(x-1)^2/6^2+y^2/3-1==0;
        E2=(x+1)^2/2^2+(y-5)^2/4^2-1==0;
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        fimplicit(E1,[-6 8 -2 10])
        hold on
        fimplicit(E2,[-6 8 -2 10])
        axis equal
        [x, y]=solve(E1,E2, 'Real', true);
        plot(x,y, 'ok', 'MarkerSize', 7, 'MarkerFaceColor','r');
        text(x(1)+1/2,y(1)+1/2,['(',num2str(round(100*double(x(1)))/100),',',num2str(round(100*double(y(1)))/100),')'])
        text(x(2)+1/2,y(2)+1/2,['(',num2str(round(100*double(x(2)))/100),',',num2str(round(100*double(y(2)))/100),')'])
        grid on
        xlabel('x', 'Interpreter', 'latex')
        ylabel('y', 'Interpreter', 'latex');
        title('10. $Elipses$', 'Interpreter', 'latex')
        clear x y E1 E2

        %L=120; CD=66; W=200;
        syms Fx Fy T Lc L W %real
        syms d positive
        E1=Fx-T*d/Lc==0;
        E2=Fy+T*sqrt(Lc^2-d^2)/Lc-W==0;
        E3=T*sqrt(Lc^2-d^2)/Lc*d-W*L==0;
        [T, Fx, Fy]=solve(E1,E2,E3,T,Fx,Fy);
        F=sqrt(Fx^2+Fy^2);
        T1=subs(T,{W,L,Lc},{200, 120, 66});
        F1=subs(F,{W,L,Lc},{200, 120, 66});
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        fplot(T1,[20 70])
        hold on
        fplot(F1,[20 70])
        d0=solve(diff(T1,d)==0);
        plot(d0,subs(T1,d0), 'ok', 'MarkerSize', 7, 'MarkerFaceColor','r');
        text(d0,subs(T1,d0)+50, ['   Min value:(',num2str(round(100*double(d0))/100),',',num2str(round(100*double(subs(T1,d0)))/100),')'])
        title('11. $Extension$ $of$ $beam$', 'Interpreter', 'latex')
        clear Fx Fy T Lc L W d E1 E2 E3 T1 F1 Fx1 Fy1 F d0
        legend('T,Newtons','F, Newtons');
        xlabel('d, inches', 'Interpreter', 'latex')
        ylabel('T//F', 'Interpreter', 'latex');
        grid on
        ylim([500 2500]);

        syms F x h mu N m g
        %m=18; h=10; mu=0.55; g=9.81;
        E1=-F*x/sqrt(x^2+h^2)+mu*N==0;
        E2=-m*g+N+F*h/sqrt(x^2+h^2)==0;
        [F,~]=solve(E1,E2,F,N);
        F1=subs(F,{m,h,mu,g},{18,10,0.55,9.81});
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        fplot(F1,[5 30]);
        x0=solve(diff(F1,x)==0);
        hold on
        plot(x0,subs(F1,x0), 'ok', 'MarkerSize', 7, 'MarkerFaceColor','r');
        text(x0,subs(F1,x0)+1, ['   Min value:(',num2str(round(100*double(x0))/100),',',num2str(round(100*double(subs(F1,x0)))/100),')'])
        xlabel('x, meters', 'Interpreter', 'latex')
        ylabel('F, Newtons', 'Interpreter', 'latex');
        title('12. $Bar$ $friction$', 'Interpreter', 'latex')
        grid on
        ylim([84 Inf])
        clear F x f mu N m g E1 E2 F N x0 F1 h

        syms x y R x0 y0
        E1=x^2+y^2==R^2;
        Y=solve(E1,y);
        Y=simplify(Y(1));
        Yf=diff(Y,x);
        Yk=simplify(y0+subs(Yf,{x,y},{x0,y0})*(x-x0));
        E2=x0^2+y0^2==R^2;
        X0=7;r=10;
        Y0=solve(subs(E2, {x0,R},{X0,r}));
        Y0=Y0(1);
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        fimplicit(subs(E1,R,r),[-11 11 -11 11]);
        hold on
        fplot(subs(Yk,{y0,x0,R},{X0,Y0,r}),[-11 11]);    
        plot([0 X0],[0 double(Y0)],'ok', 'MarkerSize', 3, 'MarkerFaceColor','b')
        plot([0, X0], [0 double(Y0)],'r--');
        xlabel('x', 'Interpreter', 'latex')
        ylabel('y', 'Interpreter', 'latex');
        title('13. $Tangent$ $to$ $circle$', 'Interpreter', 'latex')
        grid on
        axis equal
        xlim ([-11 11]);
        ylim ([-11 11]);
        clear x y x0 y0 X0 Y0 r R E1 E2 Yf Yk Y

        h=5; v=540; l=100;
        syms tet h v l t x
        %syms x positive
        P=h^2+x^2==l^2;
        x=solve(P, x);
        x=x(1);
        E=tand(tet)==(h/(x-v*t));
        Tet=solve(E,tet)+180*(1+abs(t-11.097)/(t-11.097))/2;
        V=diff(Tet,t);
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        subplot(2,1,1);
        fplot(subs(Tet, {h,l,v},{5,  100, 9}),[0 20]);
        grid on
        xlim ([0 20]);
        xlabel('t,minute','Interpreter','latex');
        ylabel('$\theta, degress$','Interpreter','latex');
        title('15. $Plane$ $over$ $radar$','Interpreter','latex');
        subplot(2,1,2);
        fplot(subs(V, {h,l,v},{5,  100, 9}),[0 20]);
        grid on
        xlim ([0 20]);
        ylim ([0 120]);
        xlabel('t,minute','Interpreter','latex');
        ylabel('${d \theta \over dt}$, $degress per minute$','Interpreter','latex');
        title('$Plane$ $over$ $radar$','Interpreter','latex');
        clear E h l P t Tet tet v V x


        syms x
        fprintf('\n16.\n   I1 = ');
        y=x^3/sqrt(1-x^2);
        i=int(y);
        disp(i);
        Y=x^2*cos(x);
        I=int(Y);
        fprintf('   I2 = ');
        disp(I);
        clear x y Y i I

        syms x
        S=cos(x)^2/(1+sin(x)^2);
        I=int(S, 0, pi);
        figure ('Units', 'Normalized', 'OuterPosition', [1/3 1/6 1/3 2/3]);
        fplot(S,[0 pi]);
        grid on
        t=text(0.5, 0.9, ['\itI = ', num2str(double(I))]);
        t.FontSize=14;
        t.Color='b';
        t.FontWeight='bold';
        t.FontName='TimesNewRoman';
        xlabel('x','Interpreter','latex');
        ylabel('y','Interpreter','latex');
        title('17. $S={cosx^2\over1+sinx^2}$','FontSize', 17,'Interpreter','latex');
        clear I S x t

        syms u v a b c
        x=a*cos(u)*sin(v);
        y=b*sin(u)*sin(v);
        z=c*cos(v);
        S=2*int(y*diff(x,u),u,pi,0);
        dz=diff(z);
        dV=S*dz;
        V=int(dV,-pi,0);
        fprintf('\n18.dV=S*dz\n');
        fprintf('   S=integral y*x from a=pi to b=0\n   S = ');
        disp(S);
        fprintf('   dz=dz(v)/dv*dv\n   dz/dv = ');
        disp(dz);
        fprintf('   V=integral dV from -pi to 0\n   dV/dv = ');
        disp(dV);
        fprintf('   V = ');
        disp(V);
        clear u v a b c x y z S dz dV V


        syms u m x A B C alp t
        u=A/sqrt(t)*exp(-x^2/(4*m*t))+B;
        fprintf('\n19.a) Left side = ');
        Left=diff(u,t);
        disp(Left);
        fprintf('      Right side = ');
        Right=diff(u,x,2)*m;
        disp(Right);
        fprintf('      Differense between them = ');
        sol=simplify(Left-Right);
        disp(sol);
        clear Left Right sol

        u=A*exp(-alp*x)*cos(alp*x-2*m*alp^2*t+B)+C;
        fprintf('   b) Left side = ');
        Left=diff(u,t);
        disp(Left);
        fprintf('      Right side = ');
        Right=diff(u,x,2)*m;
        disp(Right);
        fprintf('      Differense between them = ');
        sol=simplify(Left-Right);
        disp(sol);
        clear u m x A B C alp t Left Right sol
    end
    fprintf('\n\n\n------------------------------------');
    fprintf('\n\n\nPress 0 to exit\n\n\n');
    Chapter=input('Chapter=');
    if(Chapter==0)
        break;
    end
    close all
 end