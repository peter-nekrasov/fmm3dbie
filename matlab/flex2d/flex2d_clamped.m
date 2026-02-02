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

figure(1); clf
plot(S,rand(S.npatches,1))
hold on
plot(chnkr,'x-')
view(0,90)

V = eval_gauss(S.r);

a = 1;
b = 0.7;
c = 1/pi;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));

zk = [zk1 zk2];


% volume to volume part 

start = tic; 

A = flex2d.v2v_matgen(S,zk,1e-10);
l11 = eye(S.npts) + V.*A;

t1 = toc(start);

fprintf('%5.2e s : time to assemble v2v matrix\n',t1)



% boundary to volume part 
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
targetinfo = [];
targetinfo.r = S.r(1:2,:);
targetinfo.n = S.n(1:2,:);


start = tic; 
b2v = chunkerkernevalmat(chnkr,fkern,targetinfo);
l12 = V.*b2v;
t1 = toc(start);
fprintf('%5.2e s : time to assemble b2v matrix\n',t1)


% volume to boundary part 


targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];

start = tic; 
v2b_dir = flex2d.v2b_matgen_dir(S,zk,targinfo,1e-8);
v2b_neu = flex2d.v2b_matgen_neu(S,zk,targinfo,1e-8);
l21 = zeros(2*chnkr.npt,S.npts);
l21(1:2:end,:) = v2b_dir;
l21(2:2:end,:) = v2b_neu;

t1 = toc(start);
fprintf('%5.2e s : time to assemble v2b matrix\n',t1)


% boundary to boundary part 

fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');
kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = [];
opts.sing = 'log';

start = tic;
b2b = chunkermat(chnkr,fkern, opts);
l22 = b2b + 0.5*eye(2*chnkr.npt);
l22(2:2:end,1:2:end) = l22(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);

t1 = toc(start);
fprintf('%5.2e s : time to assemble b2b matrix\n',t1)



% form system matrix and rhs 
lhs = [l11, l12; 
    l21, l22];


rhs_vol = (4*a+2*b-c+V(:).').*sin(S.r(1,:)).*sin(S.r(2,:)) ; 
rhs_bc = zeros(chnkr.npt*2,1);
rhs_bc(1:2:end) = sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));
rhs_bc(2:2:end) = cos(chnkr.r(1,:)).*sin(chnkr.r(2,:)).*chnkr.n(1,:) ...
          + sin(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(2,:) ; 



rhs = zeros(S.npts+2*chnkr.npt,1);
rhs(1:S.npts) = rhs_vol;
rhs(S.npts+1:end) = rhs_bc;


% solve 

start = tic; 
sol = gmres(lhs,rhs,[],1e-10,100); 
t1 = toc(start);
fprintf('%5.2e s : time for dense gmres\n',t1)    

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
start1 = tic;
u = A*sol(1:S.npts)+chunkerkerneval(chnkr, ikern,...
    sol(S.npts+1:end), S.r(1:2,:));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

ref_u = (sin(S.r(1,:)).*sin(S.r(2,:))).';
err = abs(u - ref_u(:)) / max(abs(u));


figure(2); clf
scatter(S.r(1,:),S.r(2,:),8,log10(err)); 
colorbar



%%
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end

