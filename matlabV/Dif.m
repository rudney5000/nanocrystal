function U = Dif(t, V, D21, D32, m, dL, gR, g31, g21, g32, T, Tm, Am)
U  = zeros(9, 1);

Ei = Field(T, Tm, Am, t);

%Omega real and imaginary
A = Ei + gR*(V(8) + m*V(4)) + dL*(V(9) + m*V(5));
B = -dL*(V(8) + m*V(4)) + gR*(V(9) + m*V(5));

%for p11
U(1) =  2*(A.*V(8) + B.*V(9)) + 2*m*(A.*V(4) + B.*V(5)) + g21*V(2) + g31*V(3);
%for p22
U(2) =  -2*m*(A.*V(4) + B.*V(5)) - g21*V(2)+ g32*V(3);
%for p33
U(3) = -2*(A.*V(8) + B.*V(9)) - g31*V(3)- g32*V(3);

%for X21 and Y21
U(4) =  D21*V(5) + A.*(V(2) - V(1))*m + A.*V(6) + B.*V(7) - 0.5*g21*V(4);
U(5) = -D21*V(4) + B.*(V(2) - V(1))*m - A.*V(7) + B.*V(6) - 0.5*g21*V(5);
%for X32 and Y32
U(6) =  D32*V(7) - A.*V(4) - B*V(5) - m*(A.*V(6) + B.*V(7)) - 0.5*(g21 + g31+g32)*V(6);
U(7) = -D32*V(6) - B.*V(4) + A*V(5) - m*(A.*V(7) - B.*V(6)) - 0.5*(g21 + g31+g32)*V(7);
%for X31 and Y31
U(8) =  (D32 + D21)*V(9) + A.*(V(3) - V(1)) + m*(A.*V(6) - B.*V(7)) - 0.5*(g31+g32)*V(8);
U(9) = -(D32 + D21)*V(8) + B.*(V(3) - V(1)) + m*(B.*V(6) + A.*V(7)) - 0.5*(g31+g32)*V(9);
end