function [F, J] = find_F_and_J_beam(x, S, chnkr, v2v, v2b, b2v, b2b, eps, params)

mu  = x(1:S.npts);
rho = x(S.npts+1:end-1);
lam = x(end);

u = v2v*mu + b2v*rho;


p  = params.p;
kappa = params.kappa;
m  = params.m;
nu = params.nu;
c = params.c;

nV = S.npts;
nB = numel(rho);


nonlin = abs(u).^(p-1).*u;
dnonlin = p * abs(u).^(p-1);
dnonlin(u==0) = 0;
Dnonlin = spdiags(dnonlin, 0, nV, nV);


F1 = mu + m*u + eps*nu*nonlin - lam*u;

l22 = b2b + 0.5*eye(2*chnkr.npt);
l22(2:2:end,1:2:end) = l22(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);
F2 = v2b*mu + l22*rho;


F3 = sum(abs(u).^2.*S.wts(:)) - c^2;

F = [F1; F2; F3];



dF1mu  = speye(nV) + (m-lam)*v2v + eps*nu*(Dnonlin*v2v);
dF1rho = (m-lam)*b2v + eps*nu*(Dnonlin*b2v);
dF1lam = -u;

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