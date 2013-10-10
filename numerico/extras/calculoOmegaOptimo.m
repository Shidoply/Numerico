LARGO = 10;
VALOR = 30;

D = 4*VALOR*eye(LARGO);
U = -diag(VALOR*ones(1,LARGO-1),1);
L = -diag(VALOR*ones(1,LARGO-1),-1);

M = D-L-U

autovalores = eig(inv(D-L)*U);
radio = max(abs(autovalores))

wopt = 2/(1+sqrt(1-radio))
