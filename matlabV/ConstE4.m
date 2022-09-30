function ConstE4

tic

Xo11 =  1;
Xo22 =  0;
Xo33 =  0;
Xo31 =  0;
Yo31 =  0;
Xo32 =  0;
Yo32 =  0;
Xo21 =  0;
Yo21 =  0;


D21  = - 100;
D32  =   100;
D31  =  D21 + D32;

Ei   =  50;

G32  =  0.01;
G31  =  1;
G21  =  1;
GR   =  100;
DL   =  1000;
m    =  sqrt(G21);

Tm   =  50;
Tk   =  20;

R = [Xo11 Xo22 Xo33 Xo31 Yo31 Xo32 Yo32 Xo21 Yo21];

opt    = odeset('AbsTol',      1e-5, ...
                'Reltol',      1e-3);

[t, V] = ode23tb(@(t,V) ConstERPDE(V,D31,D32,D21,G31,G32,G21,GR,DL,Ei,m), [0 Tm], R, opt);

toc


ReE(:,1) = Ei + GR.*(V(:,8) + m.*V(:,4)) + DL.*(V(:,9) + m.*V(:,5));
ImE(:,1) =    - DL.*(V(:,8) + m *V(:,4)) + GR.*(V(:,9) + m.*V(:,5));
E(:,1)   = sqrt(ReE.^(2) + ImE.^(2));

N  = length(t);
dt = Tm/N;
T  = 0:dt:Tm;
U  = interp1(t,V,T,'pchip');

ReEr(:,1) = Ei + GR.*(U(:,8) + m.*U(:,4)) + DL.*(U(:,9) + m.*U(:,5));
ImEr(:,1) =    - DL.*(U(:,8) + m *U(:,4)) + GR.*(U(:,9) + m.*U(:,5));

ARE(:,1) = zeros(floor((Tm/dt)-(Tk/dt)),1);
AIE(:,1) = zeros(floor((Tm/dt)-(Tk/dt)),1);

j = 1;
for i = 1:Tm/dt
    if i >= Tk/dt
       ARE(j,1)  = ReEr(i,1);
       AIE(j,1)  = ImEr(i,1);
       j = j + 1;
    end
end

AEEA = dt.*fft(sqrt(ARE.^(2) + AIE.^(2)));
EAAE = fftshift(AEEA);
GME  = abs(EAAE);

fm = 2*pi/dt;
dv = 2*pi/(Tm - Tk);
v  = -fm/2:dv:fm/2;

toc

figure(2000)
plot(t, E(:,1))
axis([-10 Tm -5 300]);
grid on
grid minor
title('\it{The Intensity of The Field}','FontName','Arial Cyr');
xlabel('Time',                          'FontName','Arial Cyr');
ylabel('Field Intensity |E|',           'FontName','Arial Cyr');

figure(2100)
plot(v, GME(:,1))
axis([-1000 1000 -5 700]);
grid on
grid minor
title('\it{The amplitude of the Fourier image}','FontName','Arial Cyr');
xlabel('\omega',                                'FontName','Arial Cyr');
ylabel('Amplitude',                             'FontName','Arial Cyr');

figure(2200)
plot(ARE(:,1),AIE(:,1),'.')
axis([-500 500 -500 500]);
grid on
grid minor
title('\it{Time from 100 to 5000}', 'FontName','Arial Cyr');
xlabel('Re[\Omega]',               'FontName','Arial Cyr');
ylabel('Im[\Omega]',               'FontName','Arial Cyr');

toc

end