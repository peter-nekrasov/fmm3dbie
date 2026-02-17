function [F, J] = find_F_and_J_capacitor(x, S, chnkr, v2v, v2b, b2v, b2b, eps, params)

mu  = x(1:S.npts);
rho = x(S.npts+1:end-1);
lam = x(end);

u = v2v*mu + b2v*rho;

nonlin  = (1 + u).^(-2);
dnonlin = -2*(1 + u).^(-3);

delta = params.delta;
F1 = delta*mu + (1-eps)*lam*u + eps*lam*nonlin;

l22 = b2b + 0.5*eye(2*chnkr.npt);
kappa = params.kappa;
l22(2:2:end,1:2:end) = l22(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);
F2 = v2b*mu + l22*rho;

c = params.c;
F3 = sum(abs(u).^2.*S.wts(:)) - c^2;

F = [F1; F2; F3];

nV = S.npts;
nB = numel(rho);

Dnonlin = spdiags(dnonlin, 0, nV, nV);

dF1mu = delta*speye(nV) + (1-eps)*lam*v2v + eps*lam*(Dnonlin*v2v);
dF1rho = (1-eps)*lam*b2v + eps*lam*(Dnonlin*b2v);
dF1lam = (1-eps)*u + eps*nonlin;

dF2mu  = v2b;
dF2rho = l22;
dF2lam = zeros(nB,1);

Lu = (S.wts(:).*u(:)).';
dF3mu  = 2*Lu*v2v;
dF3rho = 2*Lu*b2v;
dF3lam = 0;

J = [dF1mu,  dF1rho,  dF1lam;
     dF2mu,  dF2rho,  dF2lam;
     dF3mu,  dF3rho,  dF3lam];

end