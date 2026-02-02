% solve non-linear eigen problem Ax = lam./(1+x).^2 
% using homotopy method
%
% homotopy over eps in [0,1], solve:
%                   Ax = eps*lam*x+(1-eps)*lam./(1+x).^2
%                   s.t. ||x||=1


clear; 
clc;


rng(2)

n = 100; % dimension of the problem
B = randn(n,n);
A = -B*B'-0.01*eye(n); % fake differential operator
              

 
[V,D] = eig(A);
evals = diag(D);
 
[~,idx] = mink(abs(evals),5);
idx = idx(end);
lam = evals(idx);
lam0 = lam;
x = V(:,idx);


Nstep = 50; 
opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12,...
    'MaxFunctionEvaluations',2000);

nc = 100;
c = 0;
for j=1:nc

    c = c+0.05;


    x = c*x/norm(x);
    y = [x; lam]; % initial start 
    
     
    res = zeros(Nstep+1,1);
    for i=0:Nstep
        eps = i/Nstep;
        
        fun = @(y) residual(y, A, c, eps);
        y = fsolve(fun, y, opts);
        
        x = y(1:end-1);
        lam = y(end);
        y = [x; lam];
        
        % F = 1./(x+1).^2;
        % res(i+1) = norm(A*x - lam*F);
        % fprintf('iter =%3d, eps=%8.5f, residual=%12.5e\n', ...
        %         i, eps, res(i+1));
    end
    

    F = 1./(x+1).^2;
    res(j) = norm(A*x - lam*F);
    fprintf('iter =%3d, eps=%8.5f, residual=%12.5e\n', ...
            j, eps, res(j));



end

%%
% figure
% semilogy((0:Nstep)/Nstep,res,LineWidth=2)
% grid on 
% title('$\|Ax-\lambda F(x)\|$',Interpreter='latex',FontSize=16)
% xlabel('$\epsilon$',Interpreter='latex',FontSize=16)



%%


function [r, J] = residual(y, A, c, eps)

x = y(1:end-1);
lam = y(end);

F = 1./(1+x).^2;
r = [A*x-(1-eps)*lam*x-eps*lam*F; 
    norm(x)-c]; 


n = length(x);

dr1dx = 2*eps*lam./(1+x).^3;           
Jxx = A - (1-eps)*lam*speye(n) + spdiags(dr1dx,0,n,n);
Jxl = -(1-eps)*x - eps*F;
nx = norm(x);
if nx == 0
    Jx2 = zeros(1,n);   
else
    Jx2 = (x.'/nx);     
end

J = [Jxx, Jxl;
     Jx2, 0];


end