function [L, R] = sixteenth (A, B, C)
    AB=B-A;
    BC=C-B;
    CA=A-C;
    AB=sqrt(sum(AB.^2));
    BC=sqrt(sum(BC.^2));
    CA=sqrt(sum(CA.^2));
    S=fiveteenth(A,B,C);
    R=AB*BC*CA/(4*S);
    L=2*pi*R;
end