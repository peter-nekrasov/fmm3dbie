% genpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/chunkie/')

% clear 
% clc

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
chnkr = makedatarows(chnkr,2);




figure(1); clf
plot(S,rand(S.npatches,1))
hold on
plot(chnkr,'x-')
view(0,90)

% load('../spikey_blob/geo.mat')
% S = srfr;
V = eval_gauss(S.r);

zk = 0;
nu = .1;



%% calculating curvature info



kappa = signed_curvature(chnkr);
kp = arclengthder(chnkr,kappa);
kpp = arclengthder(chnkr,kp);


chnkr = makedatarows(chnkr,2);
chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;






%% v2b

targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];


dx = chnkr.d(1,:);
dy = chnkr.d(2,:);
ds = sqrt(dx.*dx+dy.*dy);
taux = dx./ds;
tauy = dy./ds;

targinfo.du = [taux; tauy; zeros(size(taux))];

start = tic; 

[~, ipatch_id, uvs_targ, dists, flags] = get_closest_pts(S, targinfo);
v2b_free2 = bh2d.v2b_matgen_free2(S,zk,nu,targinfo,1e-2);


t1 = toc(start);
fprintf('%5.2e s : time to assemble v2b matrix\n',t1)



sum(isnan(v2b_free2(:)))


%%
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end