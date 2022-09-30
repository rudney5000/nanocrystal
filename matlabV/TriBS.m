function TriBS

Xo11 = 1;
Xo22 = 0;
Xo33 = 0;
Xo21 = 0;
Xo32 = 0;
Xo31 = 0;
Yo21 = 0;
Yo32 = 0;
Yo31 = 0;

D21  =  -30;
D32  =  30;

T    =  400;
Am   =  1;
Tm   =  2*T;

gR    =  100;
dL    =  1000;
g31  =  1;
g21  =  0.666;
g32  =  0.01;
m    =  sqrt(g21/g31);


R = [Xo11 Xo22 Xo33 Xo21 Yo21 Xo32 Yo32 Xo31 Yo31];

opt    = odeset('AbsTol',1e-7,'Reltol',1e-5);
[t, V] = ode23tb(@(t, V) Dif(t, V, D21, D32, m, dL, gR, g31, g21, g32, T, Tm, Am), [0 Tm], R, opt);

Ti(:,1) = t;
Ein(:,1) = zeros(length(Ti), 1);
for i = 1:length(t)
    Ein(i,1) = Field(T, Tm, Am, Ti(i,1));
end

ReE(:,1) = Ein(:,1) + gR*(V(:,8) + m*V(:,4)) + dL*(V(:,9) + m*V(:,5));
ImE(:,1) =            gR*(V(:,9) + m*V(:,5)) - dL*(V(:,8) + m*V(:,4));

E(:,1)   = sqrt(ReE(:,1).^(2) + ImE(:,1).^(2));

%subplot(2,2,1);
figure, plot(Ti(:,1), Ein(:,1))
axis ([0 Tm -10 500]);
%grid on
%grid minor
title('\it{External Field}', 'FontName','Arial Cyr');
xlabel('Time',               'FontName','Arial Cyr');
ylabel('E_{in}',             'FontName','Arial Cyr');

%subplot(2,2,2);
figure,plot(Ti(:,1), E(:,1))
axis ([0 Tm -10 250]);
%grid on
%grid minor
title('\it{Field Intensity }', 'FontName','Arial Cyr');
xlabel('Time',                'FontName','Arial Cyr');
ylabel('|E|',                 'FontName','Arial Cyr');

%subplot(2,2,[3,4]);
figure,plot(Ein(:,1), E(:,1))
axis([-10 500 0 500]);
%grid on
%grid minor
title('\it{Hysteresis Loop}',   'FontName','Arial Cyr');
xlabel('External Field E_{in}', 'FontName','Arial Cyr');
ylabel(' Field Intensity |E|',   'FontName','Arial Cyr');
end