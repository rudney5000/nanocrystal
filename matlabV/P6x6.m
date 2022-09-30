function D = P6x6


tic

d21  =      -400;

d32  =      400;



d31  =  d32 + d21;
g32  =  0.01;
g31  =  1;
g21  =  1;
m    =  1;
gR   =  100;
dL   =  1000;

dX   = 0.01;
XMax = 2000;

XN   = XMax/dX;
X    = 0;

N(:,:) = zeros(XN,6);
P(:,:) = zeros(XN,6);

a = 1;
b = 1;


for s = 1:XN
    [Eo,V]  = newanalitica(d31,d32,d21,g31,g32,g21,gR,dL,m,X);
    Jacob        =newmat(V,d31,d32,d21,g31,g21,gR,dL,m,X);
    LaConct    = max(real(eig(Jacob,'balance')));
    if LaConct <= 0
       N(a,1) = LaConct;
       N(a,2) = Eo;
       N(a,3) = 1- V(2) - V(1);
       N(a,4) = X;
       N(a,5) = V(1)-1+V(1)+V(2);
       N(a,6) = V(2)-1+V(1)+V(2);
       a = a + 1;
       else
       P(b,1) = LaConct;
       P(b,2) = Eo;
       P(b,3) = 1 - V(2) - V(1);
       P(b,4) = X;
       P(b,5) = V(1)-1+V(1)+V(2);
       P(b,6) = V(2)-1+V(1)+V(2);
       b = b + 1;
    end
       X = X + dX;
   
end;
  
    
N(a:XN,:) = [];
P(b:XN,:) = [];

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(52)
subplot(1,2,1);
hold on
 
plot(P(:,2),P(:,4), '.','Color', [1,0,0])
plot(N(:,2),N(:,4), '.','Color', [0,0,0])

axis([-5 1000 -5 1000]);
%title('\it{Population Difference}', 'FontName','Arial Cyr');
xlabel('|\Omega_{0}|',               'FontName','Arial Cyr');
ylabel('|\Omega|',                   'FontName','Arial Cyr');
hold off

subplot(1,2,2);
hold on
plot(P(:,1),P(:,4),'.','Color',   [1,0,0])
plot(N(:,1),N(:,4),'.','Color',   [0,0,0])
axis([-5 80 -5 1000]);
%title('\it{Lyapunov Exponent}', 'FontName','Arial Cyr');
xlabel('Max[Re[\lambda]]',       'FontName','Arial Cyr');
ylabel('|\Omega|',               'FontName','Arial Cyr');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
   



end