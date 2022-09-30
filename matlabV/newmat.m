function J = newmat(V,d31,d32,d21,g31,g21,gR,dL,m,X)
%x31x32x21 y31y32y21
J  = zeros(8,8);

J(1,1)=-g21;
J(1,2)=0;
J(1,3)=-2*V(5)*m*gR-2*V(8)*m*dL;
J(1,4)=0;
J(1,5)=-2*m*X-2*V(5)*m^2*gR-2*V(8)*m^2*dL;
J(1,6)=-2*V(8)*m*gR+2*V(5)*m*dL;
J(1,7)=0;
J(1,8)=-2*V(8)*m^2*gR+2*V(5)*m^2*dL;


J(2,1)=0;
J(2,2)=-g31;
J(2,3)=-2*X-2*V(3)*gR-2*V(6)*dL;
J(2,4)=0;
J(2,5)=-2*V(3)*m*gR-2*V(6)*m*dL;
J(2,6)=-2*V(6)*gR+2*V(3)*dL;
J(2,7)=0;
J(2,8)=-2*V(6)*m*gR+2*V(3)*m*dL;

J(3,1)=X;
J(3,2)=2*X;
J(3,3)=-(g31/2)-gR+V(4)*m*gR+V(1)*gR+2*V(2)*gR-V(7)*m*dL;
J(3,4)=m*X;
J(3,5)=-m*gR+V(4)*m^2*gR+m*V(1)*gR+2*m*V(2)*gR-V(7)*m^2*dL;
J(3,6)=-V(7)*m*gR+d31+dL-V(4)*m*dL-V(1)*dL-2*V(2)*dL;
J(3,7)=0;
J(3,8)=-V(7)*m^2*gR+m*dL-V(4)*m^2*dL-m*V(1)*dL-2*m*V(2)*dL;

J(4,1)=0;
J(4,2)=0;
J(4,3)=-m*X-V(5)*gR-V(3)*m*gR-V(8)*dL-V(6)*m*dL;
J(4,4)=-(g21/2)-g31/2;
J(4,5)=-X-V(5)*m*gR-V(3)*m^2*gR-V(8)*m*dL-V(6)*m^2*dL;
J(4,6)=-V(8)*gR-V(6)*m*gR+V(5)*dL+V(3)*m*dL;
J(4,7)=d32;
J(4,8)=-V(8)*m*gR-V(6)*m^2*gR+V(5)*m*dL+V(3)*m^2*dL;

J(5,1)=2*m*X;
J(5,2)=m*X;
J(5,3)=V(4)*gR-m*gR+2*m*V(1)*gR+m*V(2)*gR+V(7)*dL;
J(5,4)=X;
J(5,5)=-(g21/2)+V(4)*m*gR-m^2*gR+2*m^2*V(1)*gR+m^2*V(2)*gR+V(7)*m*dL;
J(5,6)=V(7)*gR-V(4)*dL+m*dL-2*m*V(1)*dL-m*V(2)*dL;
J(5,7)=0;
J(5,8)=V(7)*m*gR+d21-V(4)*m*dL+m^2*dL-2*m^2*V(1)*dL-m^2*V(2)*dL;

J(6,1)=0;
J(6,2)=0;
J(6,3)=V(7)*m*gR-d31-dL+V(4)*m*dL+V(1)*dL+2*V(2)*dL;
J(6,4)=0;
J(6,5)=V(7)*m^2*gR-m*dL+V(4)*m^2*dL+m*V(1)*dL+2*m*V(2)*dL;
J(6,6)=-(g31/2)-gR+V(4)*m*gR+V(1)*gR+2*V(2)*gR-V(7)*m*dL;
J(6,7)=m*X;
J(6,8)=-m*gR+V(4)*m^2*gR+m*V(1)*gR+2*m*V(2)*gR-V(7)*m^2*dL;

J(7,1)=0;
J(7,2)=0;
J(7,3)=V(8)*gR-V(6)*m*gR-V(5)*dL+V(3)*m*dL;
J(7,4)=-d32;
J(7,5)=V(8)*m*gR-V(6)*m^2*gR-V(5)*m*dL+V(3)*m^2*dL;
J(7,6)=-m*X-V(5)*gR+V(3)*m*gR-V(8)*dL+V(6)*m*dL;
J(7,7)=-(g21/2)-g31/2;
J(7,8)=X-V(5)*m*gR+V(3)*m^2*gR-V(8)*m*dL+V(6)*m^2*dL;

J(8,1)=0;
J(8,2)=0;
J(8,3)=-V(7)*gR+V(4)*dL-m*dL+2*m*V(1)*dL+m*V(2)*dL;
J(8,4)=0;
J(8,5)=-V(7)*m*gR-d21+V(4)*m*dL-m^2*dL+2*m^2*V(1)*dL+m^2*V(2)*dL;
J(8,6)=V(4)*gR-m*gR+2*m*V(1)*gR+m*V(2)*gR+V(7)*dL;
J(8,7)=-X;
J(8,8)=-(g21/2)+V(4)*m*gR-m^2*gR+2*m^2*V(1)*gR+m^2*V(2)*gR+V(7)*m*dL;
end