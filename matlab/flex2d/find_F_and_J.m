function [F, J] = find_F_and_J(x, S, chnkr, v2v, v2b, b2v, b2b, eps)

 
mu = x(1:S.npts);                          
rho = x(S.npts+1:end-1);                  
lam = x(end);                              
 
u = v2v*mu + b2v*rho;      

nonlin = (1 + u).^(-2);                  
dnonlin = -2*(1 + u).^(-3);                

 

 
F1 = -mu + (1-eps)*lam*u + eps*lam*nonlin;
F2 = v2b*mu - 0.5*rho + b2b*rho;
F3 = sum(abs(u).^2 .* S.wts(:))-1;
F = [F1; F2; F3];

 
nV = S.npts;
nB = chnkr.npt;

 
Dnonlin = spdiags(dnonlin, 0, nV, nV);

 
dF1mu = -speye(nV) + (1-eps)*lam*v2v + eps*lam*Dnonlin*v2v;

dF1rho = (1-eps)*lam*b2v + eps*lam*Dnonlin*b2v;

dF1lam = (1-eps)*u + eps*nonlin;

 
dF2mu = v2b;
dF2rho = -0.5*speye(nB) + b2b;
dF2lam = zeros(nB,1);

 

Lu = 2*(S.wts(:).*u(:)).';    

dF3mu = Lu*v2v;
dF3rho = Lu*b2v;
dF3lam = 0;

 
J = [dF1mu,  dF1rho,  dF1lam;
     dF2mu,  dF2rho,  dF2lam;
     dF3mu,  dF3rho,  dF3lam];

end



% function [F,J] = lap2d_find_F_and_J(x, S, chnkr, v2v, v2b, b2v, b2b)
% 
% 
% 
% mu = x(1:S.npts);
% rho = x(S.npts+1:end-1);
% lam = x(end);
% 
% 
% u = v2v*mu+b2v*rho;
% 
% F1 = -mu+(1-eps)*lam*(v2v*mu+b2v*rho)+eps*lam*Diag_u.^2;
% 
% 
% dF1mu = -eye(S.npts)+(1-eps)*lam*v2v-2*eps*lam*Diag_u.^3*v2v;
% dF1rho = (1-eps)*lam*b2v-2*eps*lam*Diag_u.^3*b2v;
% dF1lam = (1-eps)*u+eps*(1+u).^(-3);
% 
% 
% F2 = v2b*mu-.5*rho+b2b*rho;
% 
% 
% dF2mu = v2b;
% dF2rho = -.5*eye(chnkr.npt)+b2b;
% dF2lam = zeros(chnkr.npt,1);
% 
% 
% 
% Lu =2*spdiags(S.wts(:).*u(:), 0, S.npts, S.npts);
% Lu = ones(1,S.npts)*Lu;
% 
% F3 = sum(abs(u(:)).^2.*S.wts(:));
% 
% 
% dF3mu = Lu*v2v;
% dF3rho = Lu*b2v;
% dF3lam = 0;
% 
% F = [F1;F2;F3];
% 
% 
% J = [dF1mu, dF1rho, dF1lam;
%     dF2mu, dF2rho, dF2lam;
%     dF3mu, dF3rho, dF3lam];
% 
% 
% end