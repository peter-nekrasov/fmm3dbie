% genpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/chunkie/')

clear 
clc

run('/Users/yuguan/software/chunkie/startup.m')
run('/Users/yuguan/software/fmm3dbie/matlab/startup.m')
addpath '/Users/yuguan/software/fmm3dbie/src'

S = geometries.disk([],[],[4 4 4],8);

% chnkr = chunkerfunc(@(t) starfish(t,5,0,[0,0],0,1));
nch = 4*4;
cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
cparams.eps = 1e-6;
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = sort(chnkr);


rhs_vol = sin( 1.97*S.r(2,:) ) + ...
    cos(-3.52*S.r(1,:) ) + ...
    cos( 0.64*S.r(1,:) ).*sin( 4.11*S.r(2,:) ) + ...
    cos(-2.78*S.r(1,:) ).*cos(-1.23*S.r(2,:) ) + ...
    sin( 3.26*S.r(1,:) ).*sin( 0.88*S.r(2,:) ) ...
    - 1.14 ...
    + 1i*( ...
        sin(-2.69*S.r(2,:) ) + ...
        cos( 4.03*S.r(1,:) ) + ...
        cos(-1.31*S.r(1,:) ).*sin( 2.47*S.r(2,:) ) + ...
        cos( 0.92*S.r(1,:) ).*cos(-3.58*S.r(2,:) ) + ...
        sin( 2.84*S.r(1,:) ).*sin(-1.06*S.r(2,:) ) ...
        + 0.39 );


rhs_vol1 = sin(-3.08*S.r(2,:) ) + ...
    cos( 2.61*S.r(1,:) ) + ...
    cos( 4.27*S.r(1,:) ).*sin(-0.94*S.r(2,:) ) + ...
    cos(-1.86*S.r(1,:) ).*cos( 1.74*S.r(2,:) ) + ...
    sin( 0.53*S.r(1,:) ).*sin(-3.46*S.r(2,:) ) ...
    + 0.88 ...
    + 1i*( ...
        sin( 1.42*S.r(2,:) ) + ...
        cos(-4.91*S.r(1,:) ) + ...
        cos( 0.77*S.r(1,:) ).*sin( 3.19*S.r(2,:) ) + ...
        cos(-2.24*S.r(1,:) ).*cos(-1.08*S.r(2,:) ) + ...
        sin(-3.67*S.r(1,:) ).*sin( 0.66*S.r(2,:) ) ...
        - 2.06 );

rhs_bc = zeros(2*chnkr.npt,1);

rhs_bc(1:2:end) = sin(0.62*chnkr.r(2,:)) + ...
    cos(2.4*chnkr.r(1,:)) + ...
    cos(0.85*chnkr.r(1,:)).*sin(2.1*chnkr.r(2,:)) + ...
     cos(1.9*chnkr.r(1,:)).*cos(0.33*chnkr.r(2,:)) + ...
     sin(-2.6*chnkr.r(1,:)).*sin(1.45*chnkr.r(2,:))+.3 + ...
     1i*( sin(-0.23*chnkr.r(2,:)) + ...
    cos(3.64*chnkr.r(1,:)) + ...
    cos(-1.96*chnkr.r(1,:)).*sin(-0.75*chnkr.r(2,:)) + ...
     cos(0.42*chnkr.r(1,:)).*cos(-2.54*chnkr.r(2,:)) + ...
     sin(1.24*chnkr.r(1,:)).*sin(-3.26*chnkr.r(2,:))-2.4);


rhs_bc(2:2:end) = sin( 1.73*chnkr.r(2,:) ) + ...
    cos(-3.91*chnkr.r(1,:) ) + ...
    cos( 0.27*chnkr.r(1,:) ).*sin(-4.62*chnkr.r(2,:) ) + ...
    cos(-2.88*chnkr.r(1,:) ).*cos( 1.14*chnkr.r(2,:) ) + ...
    sin( 3.47*chnkr.r(1,:) ).*sin(-0.96*chnkr.r(2,:) ) ...
    - 1.83 ...
    + 1i*( ...
        sin(-2.41*chnkr.r(2,:)) + ...
        cos(4.28*chnkr.r(1,:)) + ...
        cos(-1.37*chnkr.r(1,:)).*sin(2.95*chnkr.r(2,:)) + ...
        cos(0.61*chnkr.r(1,:)).*cos(-3.72*chnkr.r(2,:)) + ...
        sin(-4.11*chnkr.r(1,:)).*sin(1.08*chnkr.r(2,:)) ...
        + 0.77 );



rhs_bc1 = zeros(2*chnkr.npt,1);

rhs_bc1(1:2:end) = sin(2.83*chnkr.r(2,:)) + ...
    cos(-4.17*chnkr.r(1,:)) + ...
    cos(1.26*chnkr.r(1,:)).*sin(-3.58*chnkr.r(2,:)) + ...
    cos(-2.41*chnkr.r(1,:)).*cos(0.97*chnkr.r(2,:)) + ...
    sin(3.92*chnkr.r(1,:)).*sin(-1.44*chnkr.r(2,:)) ...
    + 0.67 ...
    + 1i*( ...
        sin(-1.73*chnkr.r(2,:)) + ...
        cos( 4.56*chnkr.r(1,:)) + ...
        cos(-0.88*chnkr.r(1,:)).*sin( 2.14*chnkr.r(2,:)) + ...
        cos( 1.39*chnkr.r(1,:)).*cos(-3.26*chnkr.r(2,:)) + ...
        sin(-2.97*chnkr.r(1,:)).*sin( 0.81*chnkr.r(2,:)) ...
        - 1.92 );

rhs_bc1(2:2:end) = sin(-3.41*chnkr.r(2,:)) + ...
    cos(1.88*chnkr.r(1,:)) + ...
    cos(-4.05*chnkr.r(1,:)).*sin( 0.73*chnkr.r(2,:)) + ...
    cos(2.67*chnkr.r(1,:)).*cos(-1.56*chnkr.r(2,:)) + ...
    sin(0.94*chnkr.r(1,:)).*sin( 3.28*chnkr.r(2,:)) ...
    - 0.48 ...
    + 1i*( ...
        sin(2.26*chnkr.r(2,:)) + ...
        cos(-3.79*chnkr.r(1,:)) + ...
        cos(1.17*chnkr.r(1,:)).*sin(-2.64*chnkr.r(2,:)) + ...
        cos(-2.91*chnkr.r(1,:)).*cos(0.69*chnkr.r(2,:)) + ...
        sin(3.53*chnkr.r(1,:)).*sin(-1.12*chnkr.r(2,:)) ...
        + 2.31) + ...
        + chnkr.r(1,:) - 2*chnkr.r(2,:) + 1i *(0.03*chnkr.r(1,:).^2 - 2*chnkr.r(2,:));


g = [rhs_vol(:); rhs_bc(:)];
f = [rhs_vol1(:); rhs_bc1(:)];



p = 6;
% nu = 10;
nu = 1;
m = 1;



% PDE coefficient 
zk = 0;


% kernels 
v2v = flex2d.v2v_matgen(S,zk,1e-12);
v2v = real(v2v);

fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
targetinfo = [];
targetinfo.r = S.r(1:2,:);
targetinfo.n = S.n(1:2,:);
b2v = chunkerkernevalmat(chnkr,fkern,targetinfo);
b2v = real(b2v);

targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
v2b_dir = flex2d.v2b_matgen_dir(S,zk,targinfo,1e-12);
v2b_neu = flex2d.v2b_matgen_neu(S,zk,targinfo,1e-12);
v2b = zeros(2*chnkr.npt,S.npts);
v2b(1:2:end,:) = v2b_dir;
v2b(2:2:end,:) = v2b_neu;
v2b = real(v2b);


fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');
kappa = signed_curvature(chnkr);
kappa = kappa(:);
opts = [];
opts.sing = 'log';

b2b = chunkermat(chnkr,fkern, opts); 
b2b = real(b2b);
l22 = b2b + 0.5*eye(2*chnkr.npt);
l22(2:2:end,1:2:end) = l22(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);


%%
if 1==0
    a = 2;
    b = 4;
    n = 50;
    k = (0:n).';

    lam_quads = (a+b)/2+(b-a)/2*cos(pi*(n-k)/n);
    
    
    vals = zeros(n+1,1);


    for i=1:n+1
    
        l11 = eye(S.npts)+(m-lam_quads(i)^4)*v2v;
        l12 = (m-lam_quads(i)^4)*b2v;
    
        lhs = [l11, l12; 
            v2b, l22];
    
     
        sol = gmres(lhs,g,[],1e-10,100); 
      
        
        W_vol = S.wts(:);
        W_bc = zeros(2*chnkr.npt,1);
        W_bc(1:2:end) = chnkr.wts(:);
        W_bc(2:2:end) = chnkr.wts(:);
        W = [W_vol(:); W_bc(:)];
        vals(i) = 1/(sum(f.*sol.*W));
        
    
    end
    
    
    %%
    figure 
    plot(lam_quads,real(vals))
    hold on 
    plot(lam_quads,imag(vals))
    legend('real','imag')
    
    
    
    %%
    
    aa = -1;
    bb = 1;
    x = (aa+bb)/2+(bb-aa)/2*cos(pi*(n-k)/n);
    [N,X] = ndgrid(k,x);
    c2v = cos(N.*acos(X)).';
    % v2c = inv(c2v);
    cfs = c2v\vals;
    
    figure; 
    plot(log10(abs(cfs)))
    legend('|coef of cheb. poly.|')
    
    %%
    
    
    coll = zeros(n,n);
    for ii=1:(n-1)
        coll(ii,ii+1) = 1/2;
        coll(ii+1,ii) = 1/2;
    end
    
    coll(1,2) = 1;
    
    cfred = cfs(1:(end-1))/(2*cfs(end));
    
    coll(end,:)=coll(end,:)-cfred.';
    
    eis = eigs(coll,n);
    
    
    figure
    scatter(real(eis), imag(eis))
    
    eis = eis(abs(real(eis))<1);
    tol = 1E-6;
    eis = eis(abs(imag(eis))<tol);
    
    
    
    lam_quad = (eis-aa)/2*(b-a)+a;
end





%%


lam_quad = 3.203849748616112 + 0.000000000009984i;
lam = real(lam_quad^4);
l11 = eye(S.npts)+(m-lam)*v2v;
l12 = (m-lam)*b2v;
lhs = [l11, l12; 
    v2b, l22];
v = null(lhs);




%%



mu  = v(1:S.npts);
rho = v(S.npts+1:end);
u = v2v*mu + b2v*rho;




% figure; 
% clf
% scatter(S.r(1,:),S.r(2,:),8,u); 
% colorbar


%% 





c = 2.5;


params.p = p;
params.m = m;
params.nu = nu;
params.kappa = kappa;
params.c = c;




if (1==0)
    
    eps = 0.3;

    mu0  = randn(S.npts,1);
    rho0 = randn(2*chnkr.npt,1);
    lam0 = 0.8;
    x0   = [mu0; rho0; lam0];
    
    [F0, J0] = find_F_and_J_beam(x0, S, chnkr, v2v, v2b, b2v, b2b, eps, params);
    
    fprintf('||F(x0)|| = %.3e\n', norm(F0));
    
    hs = 10.^(-(2:10));
    ndir = 6;
    
    for kdir = 1:ndir
        dir = randn(size(x0));
        dir = dir / norm(dir);
    
        Jv = J0*v;
    
        errs = zeros(size(hs));
        for i = 1:numel(hs)
            h = hs(i);
            Fp = find_F_and_J_beam(x0 + h*dir, S, chnkr, v2v, v2b, b2v, b2b, eps, params);
            Fm = find_F_and_J_beam(x0 - h*dir, S, chnkr, v2v, v2b, b2v, b2b, eps, params);
    
            fd = (Fp - Fm)/(2*h);
            errs(i) = norm(fd - Jv) / max(1, norm(Jv));
        end
    
        fprintf('\nDirection %d:\n', kdir);
        fprintf('   h          relerr\n');
        for i = 1:numel(hs)
            fprintf('%8.1e   %10.3e\n', hs(i), errs(i));
        end
    
        [emin, idx] = min(errs);
        fprintf('min relerr = %.3e at h = %.1e\n', emin, hs(idx));
    end
    


end




%%

v1 = c*v/sqrt(sum(abs(u).^2.*S.wts(:)));
z = [-v1; lam]; 
z0 = z;

opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'MaxFunctionEvaluations',10000);




Nstep = 20;
for i=1:Nstep
    eps = i/Nstep;
    fun = @(y) find_F_and_J_beam(y, S, chnkr, v2v, v2b, b2v, b2b, eps, params);
    z = fsolve(fun, z, opts);
end



%%

sol = z0(1:end-1);



nx = 200;
x = linspace(-1,1,nx);
y = linspace(-1,1,nx);
[x,y] = ndgrid(x,y); 
x = x(:).';
y = y(:).';
idx = x.^2+y.^2<1;
xi = x(idx);
yi = y(idx);
targinfo = [];
targinfo.r = [xi; yi; zeros(size(xi))];
targinfo.n = zeros(size(targinfo.r));
targinfo.n(3,:) = 1;


A = flex2d.v2b_matgen_dir(S,zk,targinfo,1e-8);
ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
u = A*sol(1:S.npts)+chunkerkerneval(chnkr, ikern,...
    sol(S.npts+1:end), targinfo.r(1:2,:));
u = real(u);

U = NaN(size(x));
U(idx) = u;
U1 = reshape(U, [nx,nx]);



figure; 
tiledlayout('flow')
nexttile 
h = imagesc(U1);
set(h, 'AlphaData', ~isnan(U1));
set(gca, 'Color', 'w');  
axis image
colorbar
title('$u^{(0)} (\lambda^{(0)}\approx 105.3631)$',Interpreter='latex',FontSize=16)



sol = z(1:end-1);
u = A*sol(1:S.npts)+chunkerkerneval(chnkr, ikern,...
    sol(S.npts+1:end), targinfo.r(1:2,:));
u = real(u);


U = NaN(size(x));
U(idx) = u;
U2 = reshape(U, [nx,nx]);



nexttile 
h = imagesc(U2);
set(h, 'AlphaData', ~isnan(U2));
set(gca, 'Color', 'w');  
axis image
colorbar
title('$u^{(K)} (\lambda^{(K)}\approx 614.011)$',Interpreter='latex',FontSize=16)
