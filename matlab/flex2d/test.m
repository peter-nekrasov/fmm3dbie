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




n = 25;
k = (0:n).';
delta = 0.1;



% PDE coefficient 
a0 = -delta;
b0 = -1;
c0 = 0;
zk1 = sqrt((- b0 + sqrt(b0^2 + 4*a0*c0))/(2*a0));
zk2 = sqrt((- b0 - sqrt(b0^2 + 4*a0*c0))/(2*a0));
zk = [zk1 zk2];


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


l11 = a0*eye(S.npts);
l12 = zeros(S.npts,2*chnkr.npt);
l21 = v2b;

lhs = [l11, l12; l21, l22];

lam0 = 5;
rhs = zeros(S.npts+2*chnkr.npt,1);
rhs(1:S.npts) = lam0;



start = tic; 
sol = gmres(lhs,rhs,[],1e-10,100); 
t1 = toc(start);

%%


% ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
u = v2v*sol(1:S.npts)+b2v*sol(S.npts+1:end);
c = sqrt(sum(abs(u).^2.*S.wts(:)));



z = [sol; lam0];
z0 = z;



params = [];
params.kappa = kappa;
params.delta = delta;
params.c = c;
opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'MaxFunctionEvaluations',2000);

Nstep = 20;
for i=1:Nstep
    eps = i/Nstep;
    fun = @(y) find_F_and_J_capacitor_gdwl(y, S, chnkr, v2v, v2b, b2v, b2b, eps, params);
    z = fsolve(fun, z, opts);
end


%%

sol1 = z(1:end-1);


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




v2v_eval = real(flex2d.v2b_matgen_dir(S,zk,targinfo,1e-8));
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
b2v_eval = real(chunkerkernevalmat(chnkr,fkern,targinfo.r(1:2,:)));


u = v2v_eval*sol(1:S.npts)+b2v_eval*sol(S.npts+1:end);
u1 = v2v_eval*sol1(1:S.npts)+b2v_eval*sol1(S.npts+1:end);

%%


U = NaN(size(x));
U(idx) = u;
U = reshape(U, [nx,nx]);




U1 = NaN(size(x));
U1(idx) = u1;
U1 = reshape(U1, [nx,nx]);




figure; 
tiledlayout('flow')
nexttile 
h = imagesc(U);
set(h, 'AlphaData', ~isnan(U));
set(gca, 'Color', 'w');  
axis image
colorbar
title('$u^{(0)}$',Interpreter='latex',FontSize=16)


nexttile 
h = imagesc(U1);
set(h, 'AlphaData', ~isnan(U1));
set(gca, 'Color', 'w');  
axis image
colorbar
title('$u^{(K)}$',Interpreter='latex',FontSize=16)


nexttile 
h = imagesc(abs(U-U1));
set(h, 'AlphaData', ~isnan(U1));
set(gca, 'Color', 'w');  
axis image
colorbar
title('$|u^{(0)}-u^{(K)}|$',Interpreter='latex',FontSize=16)
