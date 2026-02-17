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
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = sort(chnkr);



zk = 0;

%v2v 
v2v = flex2d.v2v_matgen(S,zk,1e-10);



%b2v
fkern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
targetinfo = [];
targetinfo.r = S.r(1:2,:);
targetinfo.n = S.n(1:2,:);
b2v = chunkerkernevalmat(chnkr,fkern,targetinfo);


%v2b 
targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
v2b_dir = flex2d.v2b_matgen_dir(S,zk,targinfo,1e-8);
v2b_neu = flex2d.v2b_matgen_neu(S,zk,targinfo,1e-8);
l21 = zeros(2*chnkr.npt,S.npts);
l21(1:2:end,:) = v2b_dir;
l21(2:2:end,:) = v2b_neu;

%b2b
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');
kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = [];
opts.sing = 'log';

b2b = chunkermat(chnkr,fkern, opts);
l22 = b2b + 0.5*eye(2*chnkr.npt);
l22(2:2:end,1:2:end) = l22(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);




uref = sin(S.r(1,:)).*sin(S.r(2,:));
uref = uref(:);

m = 1;
lam = 0.1;
p = 6;

rhs_vol = (4+m)*sin(S.r(1,:)).*sin(S.r(2,:)) ;
rhs_vol = rhs_vol+lam*abs(sin(S.r(1,:)).*sin(S.r(2,:))).^(p-1).*sin(S.r(1,:)).*sin(S.r(2,:));

rhs_bc = zeros(chnkr.npt*2,1);
rhs_bc(1:2:end) = sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));
rhs_bc(2:2:end) = cos(chnkr.r(1,:)).*sin(chnkr.r(2,:)).*chnkr.n(1,:) ...
          + sin(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(2,:) ; 

rhs = zeros(S.npts+2*chnkr.npt,1);
rhs(1:S.npts) = rhs_vol;
rhs(S.npts+1:end) = rhs_bc;


fkern_d =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');

%%
niter = 20;
V = zeros(S.npts,1);
for iter=1:niter


    l11 = eye(S.npts)+V.*v2v;
    l12 = V.*b2v;

    lhs = [l11, l12; l21, l22];

    dens = lhs\rhs;
    mu = dens(1:S.npts);
    rho = dens(S.npts+1:end);

    u = v2v*mu + chunkerkerneval(chnkr,fkern_d,rho,S.r(1:2,:));
    u = u(:);

    rel_err = norm(u-uref)/norm(uref)
    V = m+lam*abs(u).^(p-1);
    
end





%%


