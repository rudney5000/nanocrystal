function ConstE31

tic

Xo11 =  0.3998;
Xo22 =  0.5994;
Xo33 =  0.0008;
Xo31 = -0.0000;
Yo31 = -0.0120;
Xo32 = -0.0000;
Yo32 =  0.0147;
Xo21 = -0.4891;
Yo21 =  0.0095;

D21  =  - 10;
D32  =  10;
D31  =  D21 + D32;

Ei   =  300.0000;

G32  =  0;
G31  =  1;
G21  =  0.816;
GR   =  10;
DL   =  50;
m    =  sqrt(G21);

Tm   =  1000;
Tk   =  100;

R = [Xo11 Xo22 Xo33 Xo31 Yo31 Xo32 Yo32 Xo21 Yo21];

opt    = odeset('AbsTol',      1e-8, ...
                'Reltol',      1e-6, ...
                'InitialStep', 1e-3, ...
                'MaxStep',     1e-3);

[t, V] = ode23tb(@(t,V) ECRPDE(V,D31,D32,D21,G31,G32,G21,GR,DL,Ei,m), [0 Tm], R, opt);

toc

A = Ei + gR*(V(8) + m*V(4)) + dL*(V(9) + m*V(5));
B =     -dL*(V(8) + m*V(4)) + gR*(V(9) + m*V(5));

ReE(:,1) = Ei + GR.*(V(:,8) + m.*V(:,4)) + DL.*(V(:,9) + m.*V(:,5));
ImE(:,1) =    - DL.*(V(:,8) + m *V(:,4)) + GR.*(V(:,9) + m.*V(:,5));
E(:,1)   = sqrt(ReE.^(2) + ImE.^(2));

N  = length(t);
dt = Tm/N;
T  = 0:dt:Tm;
U  = interp1(t,V,T,'pchip');

ReEr(:,1) = Ei + GR.*(U(:,8) + m.*U(:,4)) + DL.*(U(:,9) + m.*U(:,5));
ImEr(:,1) =    - DL.*(U(:,8) + m *U(:,4)) + GR.*(U(:,9) + m.*U(:,5));
Er(:,1)   = sqrt(ReEr.^(2) + ImEr.^(2));

FEEF = dt.*fft(Er);
EFFE = fftshift(FEEF);
QME  = abs(EFFE);

fm = 2*pi/dt;
dw = 2*pi/Tm;
w  = -fm/2:dw:fm/2;

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

AEEA = dt.*fft(ARE.^(2) + AIE.^(2));
EAAE = fftshift(AEEA);
GME  = abs(EAAE);

dT = (Tm - Tk)/(j-1);
fk = 2*pi/dT;
dv = 2*pi/(Tm - Tk);
v  = -fk/2:dv:fk/2;

AREk(:,1) = ThinnedVector(ARE(:,1),k);
AIEk(:,1) = ThinnedVector(AIE(:,1),k);

FREr(:,1) = ThinnedVector(ReEr(:,1),k);
FIEr(:,1) = ThinnedVector(ImEr(:,1),k);

toc

figure(201)
plot(t, E(:,1))
axis([-10 Tm -5 300]);
grid on
grid minor
title('\it{The Intensity of The Field}','FontName','Arial Cyr');
xlabel('Time',                          'FontName','Arial Cyr');
ylabel('Field Intensity |E|',           'FontName','Arial Cyr');

figure(202)
plot(w, QME(:,1))
axis([-1000 1000 -5 700]);
grid on
grid minor
title('\it{The amplitude of the Fourier image}','FontName','Arial Cyr');
xlabel('\omega',                                'FontName','Arial Cyr');
ylabel('Amplitude',                             'FontName','Arial Cyr');

figure(203)
plot(v, GME(:,1))
axis([-1000 1000 -5 700]);
grid on
grid minor
title('\it{The amplitude of the Fourier image}','FontName','Arial Cyr');
xlabel('\omega',                                'FontName','Arial Cyr');
ylabel('Amplitude',                             'FontName','Arial Cyr');

figure(204)
plot(FREr(:,1), FIEr(:,1),'.',ReE(1,1),ImE(1,1),'o')
axis([-50 50 -50 50]);
grid on
grid minor
title('\it{Phase trajectory for k = 100}','FontName','Arial Cyr');
xlabel('Re[\Omega]',                      'FontName','Arial Cyr');
ylabel('Im[\Omega]',                      'FontName','Arial Cyr');


figure(205)
plot(ReE(:,1), ImE(:,1),'.',ReE(1,1),ImE(1,1),'o')
axis([-50 50 -50 50]);
grid on
grid minor
title('\it{The full phase trajectory}', 'FontName','Arial Cyr');
xlabel('Re[\Omega]',                    'FontName','Arial Cyr');
ylabel('Im[\Omega]',                    'FontName','Arial Cyr');


figure(206)
plot(ARE(:,1),AIE(:,1),'.')
axis([-50 50 -50 50]);
grid on
grid minor
title('\it{Time from 600 to 1000}', 'FontName','Arial Cyr');
xlabel('Re[\Omega]',               'FontName','Arial Cyr');
ylabel('Im[\Omega]',               'FontName','Arial Cyr');


figure(207)
plot(AREk(:,1),AIEk(:,1),'.')
axis([-50 50 -50 50]);
grid on
grid minor
title('\it{Time from 600 to 1000}', 'FontName','Arial Cyr');
xlabel('Re[\Omega]',               'FontName','Arial Cyr');
ylabel('Im[\Omega]',               'FontName','Arial Cyr');

toc

end