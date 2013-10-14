LARGO = 70;
VALOR = 5;

D = 4*VALOR*eye(LARGO);
U = -diag(VALOR*ones(1,LARGO-1),1);
L = -diag(VALOR*ones(1,LARGO-1),-1);

%GS
M = D-L-U;
autovalores_GS = eig(inv(D-L)*U);
radio_GS = max(abs(autovalores_GS))
wopt = 2/(1+sqrt(1-radio_GS))

%Jacobi
autovalores_J = eig(inv(D)*(L+U));
radio_J = max(abs(autovalores_J))
