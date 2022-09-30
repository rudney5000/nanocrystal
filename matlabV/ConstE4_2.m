
tic

Xo11 =   0.7096;
Xo22 =   0.0000;
Xo33 =   0.2904;

Xo31 =  -0.3489;
Xo32 =   0.0000;
Xo21 =  -0.0000;
Yo31 =  -0.0008;
Yo32 =  -0.0012;
Yo21 =  -0.0024;



D21  =   -100;
D32  =   100;
D31  =  D21 + D32;

Ei   =  200;

G32  =  0.01;
G31  =  1;
G21  =  0.666;
GR   =  100;
DL   =  1000;
m    =  sqrt(G21/G31);

Tm   =  400;
Tk   =  100;

R = [Xo11 Xo22 Xo33 Xo21 Yo21 Xo32 Yo32 Xo31 Yo31 ];

opt    = odeset('AbsTol', 1e-6, ...
                'Reltol', 1e-4);

[t, V] = ode23tb(@(t,V) ConstERPDE(V,D32,D21,G31,G32,G21,GR,DL,Ei,m), [0 Tm], R, opt);

ReE(:,1) = Ei + GR.*(V(:,8) + m.*V(:,4)) + DL.*(V(:,9) + m.*V(:,5));
ImE(:,1) =    - DL.*(V(:,8) + m *V(:,4)) + GR.*(V(:,9) + m.*V(:,5));
E(:,1)   = sqrt(ReE.^(2) + ImE.^(2));

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(t);
dt = Tm/N;
T  = 0:dt:Tm;
U  = interp1(t,V,T,'pchip');

for k = 1:N
    if T(k) >= Tk
       M = k;
       break
    end
end

G(:,1) = U(M:end,4);
G(:,2) = U(M:end,5);
G(:,3) = U(M:end,6);
G(:,4) = U(M:end,7);
G(:,5) = U(M:end,8);
G(:,6) = U(M:end,9);
G(:,7) = U(M:end,3) - U(M:end,1);
G(:,8) = U(M:end,3) - U(M:end,2);
G(:,9) = U(M:end,2) - U(M:end,1);

Hm  = (2*pi)/dt;
dHz = (2*pi)/(Tm - Tk);
H   = -Hm/2:dHz:Hm/2;

REk(:,1) = Ei + GR.*(G(:,5) + m.*G(:,1)) + DL.*(G(:,6) + m.*G(:,2));
IEk(:,1) =      GR.*(G(:,6) + m.*G(:,2)) - DL.*(G(:,5) + m *G(:,1));


AEEA = (dt*N).*ifft(REk + 1i.*IEk);
EAAE = fftshift(AEEA);
GME(:,1) = H;
GME(:,2) = abs(EAAE);

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)


plot(t, E(:,1))
axis([-10 Tm -5 110]);
xlabel('t',        'FontName','Arial Cyr');
ylabel('|\Omega|', 'FontName','Arial Cyr');

figure(3)
plot(GME(:,1), GME(:,2))
axis([-1000 400 -220 6000]);
xlabel('\omega',   'FontName','Arial Cyr');
ylabel('Spectrum', 'FontName','Arial Cyr');

figure(4)
plot(REk(:,1),IEk(:,1),'.')
axis([-15 55 -20 15]);
xlabel('Re[\Omega]', 'FontName','Arial Cyr');
ylabel('Im[\Omega]', 'FontName','Arial Cyr');

toc
