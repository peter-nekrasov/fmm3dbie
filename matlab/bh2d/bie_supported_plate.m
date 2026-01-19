clear 
clc

run('/Users/yuguan/software/chunkie/startup.m')


zk = 0;
nu = 0.3;

% discretize domain

cparams = [];
cparams.eps = 1.0e-6;
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

kappa = signed_curvature(chnkr);
kp = arclengthder(chnkr,kappa);
kpp = arclengthder(chnkr,kp);

chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;



% building RHS

src = [5;5];
targ = chnkr.r(1:2,:);
[val,grad,hess] = chnk.flex2d.bhgreen(src,targ);


rhs = zeros(2*chnkr.npt, 1);
rhs(1:2:end) = val;
rhs(2:2:end) = nu*(hess(:,1)+hess(:,3))+(1-nu)*(chnkr.n(1,:).^2.'.*hess(:,1)+...
   2*(chnkr.n(1,:).*chnkr.n(2,:)).'.*hess(:,2)+chnkr.n(2,:).^2.'.*hess(:,3));


% assembling system matrix

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

sys =  M;


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
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic; 
in = chunkerinterior(chnkr,{xtarg,ytarg}); 
t1 = toc(start);

fprintf('%5.2e s : time to find points in domain\n',t1)

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_eval', nu);                              % build the kernel of evaluation          

start1 = tic;
uscat = chunkerkerneval(chnkr, ikern,sol, targets(:,in));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

utrue = chnk.flex2d.bhgreen(src,targets(:,in));
maxu = max(abs(utrue(:)));


figure(3)
clf

t = tiledlayout('flow');

nexttile
zztarg = zeros(size(xxtarg));
zztarg(in) = utrue;
% h=pcolor(xxtarg,yytarg,imag(zztarg),"FaceColor","interp");
h=pcolor(xxtarg,yytarg,real(zztarg));
shading interp
colorbar
set(h,'EdgeColor','none')
% clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{true}$','Interpreter','latex','FontSize',12)



nexttile
zztarg = zeros(size(xxtarg));
zztarg(in) = uscat;
% h=pcolor(xxtarg,yytarg,imag(zztarg),"FaceColor","interp");
h=pcolor(xxtarg,yytarg,real(zztarg));
shading interp
colorbar
set(h,'EdgeColor','none')
% clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{sol}$','Interpreter','latex','FontSize',12)



nexttile
zztarg = zeros(size(xxtarg));
zztarg(in) = log10(abs(uscat-utrue)./abs(utrue));
% h=pcolor(xxtarg,yytarg,imag(zztarg),"FaceColor","interp");
h=pcolor(xxtarg,yytarg,zztarg);
shading interp
colorbar
set(h,'EdgeColor','none')
% clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('error','Interpreter','latex','FontSize',12)


