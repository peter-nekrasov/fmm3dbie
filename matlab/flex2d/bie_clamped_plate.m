clear 
clc

run('/Users/yuguan/software/chunkie/startup.m')


a = 1;
b = 0.7;
c = 1/pi;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));


zk = [zk1 zk2];
nu = 0.3;

% discretize domain

nch = 4*4;
cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = sort(chnkr);

chnkr = makedatarows(chnkr,2);

kappa = signed_curvature(chnkr);
kp = arclengthder(chnkr,kappa);
kpp = arclengthder(chnkr,kp);

chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;




% plot geometry and data

% figure(1)
% clf
% plot(chnkr,'-x')
% hold on
% quiver(chnkr)
% axis equal


% kernel 

kern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 's');
kern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bcs');



% build RHS

src = [5;5];
srcinfo = [];
srcinfo.r = src;
targ = chnkr.r(1:2,:);


rhs = kern2(srcinfo, chnkr);


% build system matrix

fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');

kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = [];
opts.sing = 'log';

start = tic;
sys = chunkermat(chnkr,fkern, opts);
sys = sys + 0.5*eye(2*chnkr.npt);
sys(2:2:end,1:2:end) = sys(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);

% solve 

sol = gmres(sys,rhs,[],1e-12,100);


% evaluate at targets


rmin = min(chnkr); rmax = max(chnkr);
nplot = 200;
xtarg = linspace(-1.25,1.25,nplot); 
ytarg = linspace(-1.25,1.25,nplot); 
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); 
targets(2,:) = yytarg(:);
in = chunkerinterior(chnkr,{xtarg,ytarg});
targets = targets(:,in);

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
usol = chunkerkerneval(chnkr, ikern,sol, targets);

targinfo = [];
targinfo.r = targets;
utrue = kern1(srcinfo, targinfo);


% plot 

figure(1)
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
zztarg(in) = usol;
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
zztarg(in) = log10(abs(usol-utrue)./abs(utrue));
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


