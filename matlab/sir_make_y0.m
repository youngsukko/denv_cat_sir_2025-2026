function y0 = sir_make_y0(A0,A1,A2,A3,A4,R0,A5)
    mmu = 1/76.4;
    sc  = mmu/(mmu+52)*(1-1/R0);
    y0  = zeros(51*A5,1);
    y0(A5*0+(1:A5))  = A0*(1-sc);
    y0(A5*3+(1:A5))  = A1*(1-sc);
    y0(A5*6+(1:A5))  = A2*(1-sc);
    y0(A5*9+(1:A5))  = A3*(1-sc);
    y0(A5*11+(1:A5)) = A4;
    y0(A5*1+(1:A5))  = A0*sc;
    y0(A5*4+(1:A5))  = A1*sc;
    y0(A5*7+(1:A5))  = A2*sc;
    y0(A5*10+(1:A5)) = A3*sc;
end
