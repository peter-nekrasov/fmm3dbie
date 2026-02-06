
clear 
clc 


S = geometries.disk([],[],[4 4 4],8);

% chnkr = chunkerfunc(@(t) starfish(t,5,0,[0,0],0,1));
nch = 4*4;
cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = sort(chnkr);

% figure(1); clf
% plot(S,rand(S.npatches,1))
% hold on
% plot(chnkr,'x-')
% view(0,90)

 

%%

% lhs = [lhs_11, lhs_12; lhs_21, lhs_22];


%%

if 1==0

    g = zeros(S.npts+chnkr.npt,1);
    f = zeros(S.npts+chnkr.npt,1);
    
    g(1:S.npts) = sin(0.62*S.r(1,:)) + ...
        cos(2.4*S.r(1,:)) + ...
        cos(0.85*S.r(1,:)).*sin(2.1*S.r(2,:)) + ...
         cos(1.9*S.r(1,:)).*cos(0.33*S.r(2,:)) + ...
         sin(-2.6*S.r(1,:)).*sin(1.45*S.r(2,:))-1.2;
    
    
    g(S.npts+1:end) = cos(1.75*chnkr.r(2,:)) + ...
        sin(1.3*chnkr.r(1,:)).*cos(0.72*chnkr.r(2,:)) + ...
         cos(-2.2*chnkr.r(1,:)).*sin(0.58*chnkr.r(2,:)) + ...
         cos(3.1*chnkr.r(1,:)).*cos(-0.41*chnkr.r(2,:))+1;
    
    f(1:S.npts) = sin(2.83*S.r(2,:)) + ...
        cos(-4.17*S.r(1,:)) + ...
        cos(1.26*S.r(1,:)).*sin(-3.58*S.r(2,:)) + ...
        cos(-2.41*S.r(1,:)).*cos(0.97*S.r(2,:)) + ...
        sin(3.92*S.r(1,:)).*sin(-1.44*S.r(2,:)) ...
        + 0.67; 
    
    f(S.npts+1:end) = sin(-0.23*chnkr.r(2,:)) + ...
        cos(3.64*chnkr.r(1,:)) + ...
        cos(-1.96*chnkr.r(1,:)).*sin(-0.75*chnkr.r(2,:)) + ...
         cos(0.42*chnkr.r(1,:)).*cos(-2.54*chnkr.r(2,:)) + ...
         sin(1.24*chnkr.r(1,:)).*sin(-3.26*chnkr.r(2,:))-2.4; 
    
    
    
    
    
    lhs_11 = -eye(S.npts);
    lhs_12 = zeros(S.npts, chnkr.npt);
    
    
    targinfo=[];
    targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
    targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
    
    
    
    
    
    
    
    n = 20;
    k = (0:n).';
    
    
    a = 1;
    b = 10;
    zks = (a+b)/2+(b-a)/2*cos(pi*(n-k)/n);
    
    vals = zeros(n+1,1);
    for i=1:n+1
    
        lhs_21 = helm2d.v2b_dir(S,sqrt(zks(i)),targinfo,1e-8);
    
        h2d_d = kernel('h','d',sqrt(zks(i)));
        lhs_22 = -0.5*eye(chnkr.npt)+chunkermat(chnkr,h2d_d);
    
        lhs = [lhs_11, lhs_12;
            lhs_21, lhs_22];
    
    
        sol = gmres(lhs,g,[],1e-10,100); 
        W_vol = S.wts(:);
        W_bc = chnkr.wts(:);
        W = [W_vol(:); W_bc(:)];
    
        vals(i) = 1/(sum(f.*sol.*W));
    
    end
    
    
    
    
    
    
    % figure 
    % plot(zks,real(vals))
    % hold on 
    % plot(zks,imag(vals))
    % legend('real','imag')
    
    
    
    
    
    aa = -1;
    bb = 1;
    x = (aa+bb)/2+(bb-aa)/2*cos(pi*(n-k)/n);
    [N,X] = ndgrid(k,x);
    c2v = cos(N.*acos(X)).';
    % v2c = inv(c2v);
    cfs = c2v\vals;
    
    
    
    
    coll = zeros(n,n);
    for ii=1:(n-1)
        coll(ii,ii+1) = 1/2;
        coll(ii+1,ii) = 1/2;
    end
    
    coll(1,2) = 1;
    
    cfred = cfs(1:(end-1))/(2*cfs(end));
    
    coll(end,:)=coll(end,:)-cfred.';
    
    eis = eigs(coll,n);
    
    eis = eis(abs(real(eis))<1);
    tol = 1E-7;
    eis = eis(abs(imag(eis))<tol);
    
    
    zk = (eis-aa)/2*(b-a)+a
    
    
    
    
    
    
    v2v = lap2d.slp_matgen(S,1e-9);
    lhs_11 = -eye(S.npts)+zk*v2v;
    
    fkern = @(s,t) chnk.lap2d.kern(s,t,'d');
    b2v = chunkerkernevalmat(chnkr,fkern,S.r(1:2,:));
    lhs_12 = zk*b2v;
    
    targinfo=[];
    targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
    targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
    lhs_21 = lap2d.v2b_dir(S,targinfo,1e-8);
    
    b2b = chunkermat(chnkr,fkern);
    lhs_22 = -0.5*eye(chnkr.npt)+b2b;
    
    
    lhs = [lhs_11, lhs_12; 
        lhs_21, lhs_22];
    
    
    
    opts.tol = 1e-12;
    opts.maxit = 5000;
    opts.p = 50;
    
    
    [x,lambda] = eigs(lhs, 1, 0, opts);   
    % x = x/norm(x);
    
    % mu = x(1:S.npts);
    % rho = x(S.npts+1:end);
    % u = v2v*mu+b2v*rho;
end


%%

% mu = x(1:S.npts);
% rho = x(S.npts+1:end);
% u = v2v*mu+b2v*rho;

% figure 
% clf
% scatter(S.r(1,:),S.r(2,:),8,-real(u)); colorbar
% 
% 
% 
% rr = sqrt(S.r(1,:).^2+S.r(2,:).^2);
% j01 = 2.404825557695773;
% u = besselj(0, j01*rr);
% U = u(:);
% zk = j01^2;
% 
% 
% 
% figure 
% clf
% scatter(S.r(1,:),S.r(2,:),8,real(U)); colorbar
% 
% 
% 
% 
% u = real(u(:));
% U = real(U(:));
% corr(u,U)

%%

zk = 5.783185962946785;
v2v = lap2d.slp_matgen(S,1e-9);
lhs_11 = -eye(S.npts)+zk*v2v;

fkern = @(s,t) chnk.lap2d.kern(s,t,'d');
b2v = chunkerkernevalmat(chnkr,fkern,S.r(1:2,:));
lhs_12 = zk*b2v;

targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
v2b = lap2d.v2b_dir(S,targinfo,1e-8);
lhs_21 = v2b;

b2b = chunkermat(chnkr,fkern);
lhs_22 = -0.5*eye(chnkr.npt)+b2b;


lhs = [lhs_11, lhs_12; 
    lhs_21, lhs_22];


opts = [];
opts.tol = 1e-12;
opts.maxit = 5000;
opts.p = 50;


[x,lambda] = eigs(lhs, 1, 0, opts);   
mu = x(1:S.npts);
rho = x(S.npts+1:end);

eps = 0.000;
lam = zk;

u = v2v*mu + b2v*rho;   

x = x/sqrt(sum(abs(u).^2 .* S.wts(:)));
mu = x(1:S.npts);
rho = x(S.npts+1:end);
u = v2v*mu + b2v*rho; 


nonlin = (1 + u).^(-2);                  
dnonlin = -2*(1 + u).^(-3);                

F1 = -mu + (1-eps)*lam*u;
F2 = v2b*mu - 0.5*rho + b2b*rho;
F3 = sum(abs(u).^2 .* S.wts(:))-1;





%%

y = [-x;lam];

Nstep = 100; 
opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'MaxFunctionEvaluations',2000);

for i=1:Nstep
    eps = i/Nstep;
    fun = @(x) find_F_and_J(x, S, chnkr, v2v, v2b, b2v, b2b, eps);
    y = fsolve(fun, y, opts);
end


%%

x = y(1:end-1);
mu = x(1:S.npts);
rho = x(S.npts+1:end);
u = v2v*mu + b2v*rho;   


figure 
scatter(S.r(1,:),S.r(2,:),8,-real(u)); colorbar