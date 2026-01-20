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

v2b_supp2 = bh2d.v2b_matgen_supp2(S,zk,nu,targinfo,1e-8);
v2b_free2 = bh2d.v2b_matgen_free2(S,zk,nu,targinfo,1e-8);
v2b = zeros(2*chnkr.npt,S.npts);
v2b(1:2:end,:) = v2b_supp2;
v2b(2:2:end,:) = v2b_free2;

t1 = toc(start);
fprintf('%5.2e s : time to assemble v2b matrix\n',t1)



%% b2b
fkern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate', nu);        % build the desired kernel
double = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';



start = tic;
sysmat1 = chunkermat(chnkr,fkern1, opts);
D = chunkermat(chnkr, double, opts);
H = chunkermat(chnkr, hilbert, opts2);  

sysmat = zeros(2*chnkr.npt);
sysmat(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) - sysmat1(3:4:end,1:2:end)*H  + 2*((1+nu)/2)^2*D*D;
sysmat(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) - sysmat1(4:4:end,1:2:end)*H;
sysmat(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
sysmat(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);

D = -[-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];  % jump matrix 
D = kron(eye(chnkr.npt), D);

b2b =  D + sysmat;

t1 = toc(start);
fprintf('%5.2e s : time to assemble b2b matrix\n',t1)



%% b2v


fkern =  @(s,t) chnk.flex2d.kern(0, s, t, 'free_plate_eval', nu);
targetinfo = [];
targetinfo.r = S.r(1:2,:);
targetinfo.n = S.n(1:2,:);


b2v = zeros([S.npts, chnkr.npt*2]);



start = tic; 
out = chunkerkernevalmat(chnkr,fkern,targetinfo);
% b2v = V.*b2v;



K1 = out(:,1:3:end);
K1H = out(:,2:3:end);
K2 = out(:,3:3:end);


b2v(:,1:2:end) = K1-K1H*H;
b2v(:,2:2:end) = K2;
b2v = V.*b2v;



t1 = toc(start);
fprintf('%5.2e s : time to assemble b2v matrix\n',t1)


%% build system matrix

sys = [v2v, b2v; 
    v2b, b2b];

% sys = @(x) [lhs_11(x(1:S.npts)) + b2v*x(S.npts+1:end); v2b*x(1:S.npts) + b2b*x(S.npts+1:end)];


%% form right hand side


rhs_vol = (4+V(:).').*sin(S.r(1,:)).*sin(S.r(2,:)); 


sinx = sin(chnkr.r(1,:));
siny = sin(chnkr.r(2,:));
cosx = cos(chnkr.r(1,:));
cosy = cos(chnkr.r(2,:));



rhs_bc = zeros(chnkr.npt*2,1);
rhs_bc(1:2:end) = -(1+nu)*siny.*siny + ...
    2*(1-nu)*cosx.*cosy.*chnkr.n(1,:).*chnkr.n(2,:);


nx = chnkr.n(1,:).';
ny = chnkr.n(2,:).';


gsxx = -sinx.*siny; 
gsxx = gsxx.';
gsxy = cosx.*siny; 
gsxy = gsxy.';
gsyy = -sinx.*siny; 
gsyy = gsyy.';
gsxxx = -cosx.*siny; 
gsxxx = gsxxx.';
gsxxy = -sinx.*cosy; 
gsxxy = gsxxy.';
gsxyy = -cosx.*siny; 
gsxyy = gsxyy.';
gsyyy = -sinx.*cosy; 
gsyyy = gsyyy.';
kappa = kappa(:);
taux = taux.';
tauy = tauy.';

rhs_bc(2:2:end) = ...
( gsxxx.*(nx.*nx.*nx) + gsxxy.*(3*nx.*nx.*ny) + gsxyy.*(3*nx.*ny.*ny) + ...
gsyyy.*(ny.*ny.*ny) ) + (2-nu) .* ( gsxxx.*(taux.*taux.*nx) + ...
           gsxxy.*(taux.*taux.*ny + 2*taux.*tauy.*nx) + ...
           gsxyy.*(2*taux.*tauy.*ny + tauy.*tauy.*nx) + ...
           gsyyy.*(tauy.*tauy.*ny) ) + ...
(1-nu) .* kappa .* ( ...
    ( gsxx.*(taux.*taux) + gsxy.*(2*taux.*tauy) + gsyy.*(tauy.*tauy) ) - ...
    ( gsxx.*(nx.*nx)     + gsxy.*(2*nx.*ny)     + gsyy.*(ny.*ny) ) );



rhs = zeros(S.npts+2*chnkr.npt,1);
rhs(1:S.npts) = rhs_vol;
rhs(S.npts+1:end) = rhs_bc;


%% solve 

start = tic; 
sol = gmres(sys,rhs,[],1e-10,200); 
t1 = toc(start);
fprintf('%5.2e s : time for dense gmres\n',t1)    


% %% postprocess 
% 
% 
% ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu);
% start1 = tic;
% u = A*(sol(1:S.npts))+chunkerkerneval(chnkr, ikern,...
%     sol(S.npts+1:end), S.r(1:2,:));
% t2 = toc(start1);
% fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)
% 
% 
% ref_u = sin(S.r(1,:)).*sin(S.r(2,:));
% err = abs(u - ref_u(:)) / max(abs(u));
% 
% 
% 
% %% ploting 
% 
% figure(2); clf
% scatter(S.r(1,:),S.r(2,:),8,log10(err)); 
% colorbar


%%
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end