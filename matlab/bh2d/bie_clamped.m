clear 
clc

run('/Users/yuguan/software/chunkie/startup.m')


zk = 0;

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
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);


fprintf('%5.2e s : time to build geo\n',t1)


% plot geometry and data

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
axis equal



% building RHS
src = [5;5];
targ = chnkr.r(1:2,:);
[val,grad] = chnk.flex2d.bhgreen(src,targ);


nx = chnkr.n(1,:); 
ny = chnkr.n(2,:);
normalderiv = grad(:, 1).*(nx.')+ grad(:, 2).*(ny.'); 

rhs = zeros(2*chnkr.npt, 1);
rhs(1:2:end) = val;
rhs(2:2:end) = normalderiv;



% assembling system matrix

fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');

kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = [];
opts.sing = 'log';

start = tic;
sys = chunkermat(chnkr,fkern, opts);
sys = sys + 0.5*eye(2*chnkr.npt);
sys(2:2:end,1:2:end) = sys(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);

t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

% Solving linear system

start = tic; sol = gmres(sys,rhs,[],1e-12,100); t1 = toc(start);
fprintf('%5.2e s : time for dense gmres\n',t1)    



% evaluate at targets and plot

rmin = min(chnkr); rmax = max(chnkr);
% xl = rmax(1)-rmin(1);
% yl = rmax(2)-rmin(2);
nplot = 200;
% xtarg = linspace(rmin(1)-xl/2,rmax(1)+xl/2,nplot); 
% ytarg = linspace(rmin(2)-yl/2,rmax(2)+yl/2,nplot);
% [xxtarg,yytarg] = meshgrid(xtarg,ytarg);
xtarg = linspace(-1.25,1.25,nplot); 
ytarg = linspace(-1.25,1.25,nplot); 
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic; in = chunkerinterior(chnkr,{xtarg,ytarg}); t1 = toc(start);

fprintf('%5.2e s : time to find points in domain\n',t1)

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');                              % build the kernel of evaluation          

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


