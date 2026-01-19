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

V = eval_gauss(S.r);

zk = 0;
nu = .1;



%% calculating curvature info

kappa = signed_curvature(chnkr);
kp = arclengthder(chnkr,kappa);
kpp = arclengthder(chnkr,kp);

chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;




%% v2v

start = tic; 
A = bh2d.v2v_matgen(S,zk,1e-8);
v2v = eye(S.npts) + V.*A;
t1 = toc(start);
fprintf('%5.2e s : time to assemble v2v matrix\n',t1)


%% b2v


fkern =  @(s,t) chnk.flex2d.kern(0, s, t, 'supported_plate_eval', nu);
targetinfo = [];
targetinfo.r = S.r(1:2,:);
targetinfo.n = S.n(1:2,:);

start = tic; 
b2v = chunkerkernevalmat(chnkr,fkern,targetinfo);
b2v = V.*b2v;
t1 = toc(start);
fprintf('%5.2e s : time to assemble b2v matrix\n',t1)


%% v2b

targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];

start = tic; 

v2b_dir = bh2d.v2b_matgen_dir(S,zk,targinfo,1e-8);
v2b_supp2 = bh2d.v2b_matgen_supp2(S,zk,nu,targinfo,1e-8);
v2b = zeros(2*chnkr.npt,S.npts);
v2b(1:2:end,:) = v2b_dir;
v2b(2:2:end,:) = v2b_supp2;

t1 = toc(start);
fprintf('%5.2e s : time to assemble v2b matrix\n',t1)



%% b2b
fkern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_log',nu);           % build the desired kernel
fkern2 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_smooth',nu);           % build the desired kernel


kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.quad = 'native';
opts2.sing = 'smooth';


start = tic;

M = chunkermat(chnkr,fkern1, opts);
M2 = chunkermat(chnkr,fkern2, opts2);

c0 = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + M2 - c0.*kappa(:).^2.*eye(chnkr.npt) - max(zk.^2)/2*eye(chnkr.npt); % extra term shows up for the general problem
M = M + 0.5*eye(2*chnkr.npt);


b2b = M;

t1 = toc(start);
fprintf('%5.2e s : time to assemble b2b matrix\n',t1)


%% form system matrix

sys = [v2v, b2v; 
    v2b, b2b];

% sys = @(x) [lhs_11(x(1:S.npts)) + b2v*x(S.npts+1:end); v2b*x(1:S.npts) + b2b*x(S.npts+1:end)];


%% form right hand side


rhs_vol = (4+V(:).').*sin(S.r(1,:)).*sin(S.r(2,:)) ; 
rhs_bc = zeros(chnkr.npt*2,1);
rhs_bc(1:2:end) = sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));
% rhs_bc(2:2:end) = cos(chnkr.r(1,:)).*sin(chnkr.r(2,:)).*chnkr.n(1,:) ...
%           + sin(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(2,:) ; % sin(xs(1,:))
rhs_bc(2:2:end) = -(1+nu)*sin(chnkr.r(1,:)).*sin(chnkr.r(2,:)) + ...
    2*(1-nu)*cos(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(1,:).*chnkr.n(2,:);



rhs = zeros(S.npts+2*chnkr.npt,1);
rhs(1:S.npts) = rhs_vol;
rhs(S.npts+1:end) = rhs_bc;


%% solve 

start = tic; 
sol = gmres(sys,rhs,[],1e-10,100); 
t1 = toc(start);
fprintf('%5.2e s : time for dense gmres\n',t1)    


%% postprocess 


ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_eval', nu);
start1 = tic;
u = A*(sol(1:S.npts))+chunkerkerneval(chnkr, ikern,...
    sol(S.npts+1:end), S.r(1:2,:));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)


ref_u = sin(S.r(1,:)).*sin(S.r(2,:));
err = abs(u - ref_u(:)) / max(abs(u));



%% ploting 

figure(2); clf
scatter(S.r(1,:),S.r(2,:),8,log10(err)); 
colorbar


%%
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end