function D = Test_2

tic

D21  = -100;
D32  =  100;
D31  =  D32 + D21;


Ei   =  350;

G32  =  0.01;
G31  =  1;
G21  =  0.666;
GR   =  100;
DL   =  1000;
m    =  sqrt(G21/G31);

dX   = 0.05;
XMax = 600;
XN   = XMax/dX;
X    = 0;

W(:,:) = zeros(XN,2);

for k = 1:XN
    Ea = newanalitica(D31,D32,D21,G31,G32,G21,GR,DL,m,X);
    W(k,1) = Ea;
    W(k,2) = X;
    X = X + dX;
end

i = 1;
U(:,:) = zeros(4,2);
P(:,1) = zeros(1,1);

for k = 4:XN
    if (W(k-1,1)>=Ei) && (W(k-2,1)<=Ei)
        U(:,:) = W(k-3:k,:);
        P(i,1) = interp1(U(:,1),U(:,2),Ei,'pchip');
        i = i + 1;
    elseif (W(k-1,1)<=Ei) && (W(k-2,1)>=Ei)
        U(:,:) = W(k-3:k,:);
        P(i,1) = interp1(U(:,1),U(:,2),Ei,'pchip');
        i = i + 1;
    end
end

D(1:i-1,1) = Ei;
D(:,2) = P(:,1);


figure(1)
plot(W(:,1),W(:,2),'-',D(:,1),D(:,2),'o') %'Color',[0,0.4470,0.7410],'-o',,'MarkerEdgeColor','b'
axis([-5 100 -5 100]);
grid on
grid minor
%title('\it{Population Difference}', 'FontName','Arial Cyr');
xlabel('|\Omega_{0}|',               'FontName','Arial Cyr');
ylabel('|\Omega|',                   'FontName','Arial Cyr');

toc

end