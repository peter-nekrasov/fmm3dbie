clear 
clc

run('/Users/yuguan/software/chunkie/startup.m')


zk = 0;
nu = 0.3;

% discretize domain

cparams = [];
cparams.eps = 1.0e-12;
cparams.nover = 0;
cparams.maxchunklen = .25; % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.25;

start = tic; 
chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
chnkr = makedatarows(chnkr,2);
t1 = toc(start);


fprintf('%5.2e s : time to build geo\n',t1)


% plot geometry and data

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
axis equal


% calculating curvature info

% kappa = signed_curvature(chnkr);
% kp = arclengthder(chnkr,kappa);
% kpp = arclengthder(chnkr,kp);
% 
% chnkr.data(1,:,:) = kp;
% chnkr.data(2,:,:) = kpp;
% 
% kappa = kappa(:);



% building RHS

src = [5;5];
targ = chnkr.r(1:2,:);
[val,der,der2,der3] = chnk.flex2d.bhgreen(src,targ);




% dx = chnkr.d(1,:).';
% dy = chnkr.d(2,:).';
% nx = chnkr.n(1,:).';
% ny = chnkr.n(2,:).';
% ds = sqrt(dx.*dx+dy.*dy);
% taux = dx./ds;
% tauy = dy./ds;
% gsxx = der2(:,1);
% gsxy = der2(:,2);
% gsyy = der2(:,3);
% gsxxx = der3(:,1);
% gsxxy = der3(:,2);
% gsxyy = der3(:,3);
% gsyyy = der3(:,4);
% 
% 
% rhs = zeros(2*chnkr.npt, 1);
% rhs(1:2:end) = nu.*(gsxx + gsyy) + ...
%     (1 - nu)*(nx.*nx.*gsxx + 2*nx.*ny.*gsxy + ny.*ny.*gsyy);
% 
% rhs(2:2:end) = ...
%     ( gsxxx.*(nx.*nx.*nx) + gsxxy.*(3*nx.*nx.*ny) + gsxyy.*(3*nx.*ny.*ny) + gsyyy.*(ny.*ny.*ny) ) + ...
%     (2-nu) .* ( gsxxx.*(taux.*taux.*nx) + ...
%                gsxxy.*(taux.*taux.*ny + 2*taux.*tauy.*nx) + ...
%                gsxyy.*(2*taux.*tauy.*ny + tauy.*tauy.*nx) + ...
%                gsyyy.*(tauy.*tauy.*ny) ) + ...
%     (1-nu) .* kappa .* ( ...
%         ( gsxx.*(taux.*taux) + gsxy.*(2*taux.*tauy) + gsyy.*(tauy.*tauy) ) - ...
%         ( gsxx.*(nx.*nx)     + gsxy.*(2*nx.*ny)     + gsyy.*(ny.*ny) ) );


srcinfo = []; 
srcinfo.r = [5;5]; 
kern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_bcs',nu);
rhs = kern2(srcinfo,chnkr);

% assembling system matrix

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

sys =  D + sysmat;

v = zeros(2*chnkr.npt,3);
w = zeros(2*chnkr.npt,3);

v(1:2:end,1) = cos(1.7*chnkr.r(1,:)) + exp(-chnkr.r(2,:)+0.3);
v(2:2:end,1) = -.3*exp(-2*chnkr.r(2,:));
v(1:2:end,2) = cos(chnkr.r(1,:).^2+0.1) + chnkr.r(1,:).*chnkr.r(2,:);
v(2:2:end,2) = sin(1.3*chnkr.r(2,:)) + exp(-chnkr.r(1,:)-0.3);
v(1:2:end,3) = cos(2*chnkr.r(2,:)+0.4) + exp(-chnkr.r(1,:));
v(2:2:end,3) = sin(chnkr.r(2,:).^2) + exp(-chnkr.r(1,:).^2);

w(1:2:end,3) = exp(1.1*chnkr.r(1,:)) + cos(-chnkr.r(2,:)+0.3);
w(2:2:end,3) = sin(-2.2*chnkr.r(2,:));
w(1:2:end,2) = exp(0.1*chnkr.r(1,:).^2+0.3) + sin(chnkr.r(1,:)).*chnkr.r(2,:);
w(2:2:end,2) = sin(1.3*chnkr.r(2,:)) - exp(-chnkr.r(1,:)-0.3);
w(1:2:end,1) = sin(2.1*chnkr.r(2,:)+0.7) + exp(-chnkr.r(1,:).^4);
w(2:2:end,1) = -exp(chnkr.r(2,:).^2) + cos(-chnkr.r(1,:).^2);

sys = sys + v*w.';

t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

% Solving linear system

start = tic; 
sol = gmres(sys,rhs,[],1e-12,100); 
t1 = toc(start);
fprintf('%5.2e s : time for dense gmres\n',t1)    



% evaluate at targets and plot

rmin = min(chnkr); rmax = max(chnkr);
nplot = 200;

xtarg = linspace(-1.25,1.25,nplot); 
ytarg = linspace(-1.25,1.25,nplot); 
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); 
targets(2,:) = yytarg(:);

start = tic; 
in = chunkerinterior(chnkr,{xtarg,ytarg}); 
targets = targets(:,in);
t1 = toc(start);

fprintf('%5.2e s : time to find points in domain\n',t1)

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu);                              % build the kernel of evaluation          

start1 = tic;

dens_comb = zeros(3*chnkr.npt,1);
dens_comb(1:3:end) = sol(1:2:end);
dens_comb(2:3:end) = -H*sol(1:2:end);
dens_comb(3:3:end) = sol(2:2:end);


Dsol = chunkerkerneval(chnkr, ikern, dens_comb, targets);
t2 = toc(start1);


fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)


utrue = chnk.flex2d.bhgreen(src,targets);
maxu = max(abs(utrue(:)));


% postprocessing 
Amat = zeros(3);
Amat(:,1) = 1;
Amat(:,2) = targets(1,1:3);
Amat(:,3) = targets(2,1:3);

cs = Amat \ (utrue(1:3) - Dsol(1:3));

uscat = Dsol + cs(1)*ones(numel(Dsol),1) + cs(2)*targets(1,:).' + cs(3)*targets(2,:).';

% ploting

figure(3)
clf

t = tiledlayout('flow');

nexttile
zztarg = zeros(size(xxtarg));
zztarg(in) = utrue;
h=pcolor(xxtarg,yytarg,real(zztarg));
shading interp
colorbar
set(h,'EdgeColor','none')
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{true}$','Interpreter','latex','FontSize',12)



nexttile
zztarg = zeros(size(xxtarg));
zztarg(in) = uscat;
h=pcolor(xxtarg,yytarg,real(zztarg));
shading interp
colorbar
set(h,'EdgeColor','none')
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{sol}$','Interpreter','latex','FontSize',12)



nexttile
zztarg = zeros(size(xxtarg));
zztarg(in) = log10(abs(uscat-utrue)./abs(utrue));
h=pcolor(xxtarg,yytarg,zztarg);
shading interp
colorbar
set(h,'EdgeColor','none')
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('error','Interpreter','latex','FontSize',12)


