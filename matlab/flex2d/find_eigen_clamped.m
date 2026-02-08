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



n = 25;
k = (0:n).';
delta = 0.1;

a = 2;
b = 2.5;
lam_quads = (a+b)/2+(b-a)/2*cos(pi*(n-k)/n);


vals = zeros(n+1,1);
l2_dens = zeros(n+1,1);


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

for i=1:n+1

    lam = -lam_quads(i)^4;
    l11 = -delta*eye(S.npts)-lam*v2v;
    l12 = -lam*b2v;

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

%%

lam_quad = real((eis-aa)/2*(b-a)+a);
lam = -lam_quad^4;
l11 = -delta*eye(S.npts)-lam*v2v;
l12 = -lam*b2v;

lhs = [l11, l12; 
    v2b, l22];
v = null(lhs);

%%
mu  = v(1:S.npts);
rho = v(S.npts+1:end);
u = v2v*mu + b2v*rho;
c = 0.1;
v = c*v/sqrt(sum(abs(u).^2.*S.wts(:)));
y = [-v; lam];

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

Nstep = 100;
for i=1:Nstep
    eps = i/Nstep;
    fun = @(y) find_F_and_J(y, S, chnkr, v2v, v2b, b2v, b2b, eps, params);
    y = fsolve(fun, y, opts);
end


%%

sol = y(1:end-1);
mu  = sol(1:S.npts);



rho = sol(S.npts+1:end);
u = v2v*mu + b2v*rho;

figure; clf
scatter(S.r(1,:),S.r(2,:),8,u); 
colorbar


%%




% rs = linspace(0,1,50);
% ths = linspace(0,2*pi,50);
% [rs, ths] = meshgrid(rs, ths);
% xs = rs.*cos(ths);
% ys = rs.*sin(ths);
% 
% 
% mu = y(1:S.npts);
% rho = y(S.npts+1:end-1);
% lam = y(end);
% fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
% targetinfo = [];
% targetinfo.r = [xs(:).';ys(:).'];
% targetinfo.n = zeros(2,50^2);
% u2 = real(chunkerkernevalmat(chnkr,fkern,targetinfo))*rho;
% 
% targetinfo.r = [xs(:).';ys(:).';zeros(1,50^2)];
% targetinfo.n = [zeros(2,50^2);-ones(1,50^2)];
% u1 =  real(flex2d.v2b_matgen_dir(S,zk,targetinfo,1e-12))*mu;
% u3 = u1+u2;
% 
% 
% 
% Nr  = size(xs,1);
% Nth = size(xs,2);
% 
% u3 = reshape(u3, Nr, Nth);
% 
% figure
% surf(xs, ys, u3, 'EdgeColor','none');
% axis equal tight
% view(2)
% colorbar