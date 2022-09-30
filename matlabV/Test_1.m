function V = Test_1

tic

D21  = -100;
D32  =  100;
D31  =  D32 + D21;

G32  =  0.01;
G31  =  1;
G21  =  0.666;
GR   =  100;
DL   =  1000;
m    =  sqrt(G21/G31);

X = 0.4161;

[Eo,V] = newanalitica9x9(D31,D32,D21,G31,G32,G21,GR,DL,m,X);

Ei = Eo;

toc

end