function [F, J] = find_F_and_J_capacitor_gdwl(x, S, chnkr, v2v, v2b, b2v, b2b, eps, params)

mu  = x(1:S.npts);
rho = x(S.npts+1:end-1);
lam = x(end);

u = v2v*mu + b2v*rho;

nonlin  = (1 + u).^(-2);
dnonlin = -2*(1 + u).^(-3);

delta = params.delta;

% Updated F1
nV = S.npts;
F1 = delta*mu - (1-eps)*lam*ones(nV,1) - eps*lam*nonlin;

l22 = b2b + 0.5*eye(2*chnkr.npt);
kappa = params.kappa;
l22(2:2:end,1:2:end) = l22(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);
F2 = v2b*mu + l22*rho;

c = params.c;
F3 = sum(abs(u).^2 .* S.wts(:)) - c^2;

F = [F1; F2; F3];

nB = numel(rho);

Dnonlin = spdiags(dnonlin, 0, nV, nV);

% Jacobians for updated F1
dF1mu  = delta*speye(nV) - eps*lam*(Dnonlin*v2v);
dF1rho =              - eps*lam*(Dnonlin*b2v);
dF1lam = -(1-eps)*ones(nV,1) - eps*nonlin;

% Jacobians for F2
dF2mu  = v2b;
dF2rho = l22;
dF2lam = zeros(nB,1);

% Jacobians for F3 (correct for complex u with real x)
Lu = (S.wts(:) .* conj(u(:))).';
dF3mu  = 2*real(Lu*v2v);
dF3rho = 2*real(Lu*b2v);
dF3lam = 0;

J = [dF1mu,  dF1rho,  dF1lam;
     dF2mu,  dF2rho,  dF2lam;
     dF3mu,  dF3rho,  dF3lam];

end