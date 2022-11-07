function th = twelwth (A, B, C)
    BA=A-B;
    BC=C-B;
    th=acosd(sum(BA.*BC)/(2*sqrt(sum(BA.^2))*sqrt(sum(BC.^2))));
end