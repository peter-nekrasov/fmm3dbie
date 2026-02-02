% genpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/chunkie/')
clear 
clc 


run('/Users/yuguan/software/chunkie/startup.m')
run('/Users/yuguan/Dropbox/fmm3dbie/matlab/startup.m')
addpath '/Users/yuguan/Dropbox/fmm3dbie/src'

S = geometries.disk([],[],[4 4 4],8);

% chnkr = chunkerfunc(@(t) starfish(t,5,0,[0,0],0,1));
nch = 4*4;
cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = sort(chnkr);








v2v = lap2d.slp_matgen(S,1e-9);

fkern = @(s,t) chnk.lap2d.kern(s,t,'d');
b2v = chunkerkernevalmat(chnkr,fkern,S.r(1:2,:));

targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
v2b = lap2d.v2b_dir(S,targinfo,1e-8);

l2d_d = kernel('l','d');
b2b = chunkermat(chnkr,l2d_d); 

fkern_d = kernel('l','d');




uref = sin(S.r(1,:)).*sin(S.r(2,:));
uref = uref(:);

a = 100;
f = sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));
g = -2*sin(S.r(1,:)).*sin(S.r(2,:)) ;
g = g+a*abs(sin(S.r(1,:)).*sin(S.r(2,:))).^2.*sin(S.r(1,:)).*sin(S.r(2,:));
rhs = [g(:); f(:)];


%%
niter = 20;
V = zeros(S.npts,1);
for iter=1:niter


    l11 = -eye(S.npts)+V.*v2v;
    l12 = V.*b2v;
    l21 = v2b;
    l22 = -0.5*eye(chnkr.npt)+b2b;

    lhs = [l11, l12; l21, l22];

    dens = lhs\rhs;
    mu = dens(1:S.npts);
    rho = dens(S.npts+1:end);

    u = v2v*mu + chunkerkerneval(chnkr,fkern_d,rho,S.r(1:2,:));
    u = u(:);

    rel_err = norm(u-uref)/norm(uref)
    V = a*abs(u).^2;
    
end