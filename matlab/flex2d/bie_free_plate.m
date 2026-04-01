a = 1;
b = 0.7;
c = 1/pi;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));


zk = [zk1 zk2];
nu = 0.3;

% discretize domain

cparams = []; cparams.maxchunklen = 4/max(abs(zk));
chnkr = chunkerfunc(@(t) ellipse(t),cparams);
chnkr = sort(chnkr);



% plot geometry and data

% figure(1)
% clf
% plot(chnkr,'-x')
% hold on
% quiver(chnkr)
% axis equal


% kernel 

kern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 's');
kern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_bcs', nu);



% build RHS

src = [5;5];
srcinfo = [];
srcinfo.r = src;
targ = chnkr.r(1:2,:);


rhs = kern2(srcinfo, chnkr);


% build system matrix

fkern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate', nu);        % build the desired kernel
double = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';

% building system matrix

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

t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

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

dens_comb = zeros(3*chnkr.npt,1);
dens_comb(1:3:end) = sol(1:2:end);
dens_comb(2:3:end) = -H*sol(1:2:end);
dens_comb(3:3:end) = sol(2:2:end);

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu);
usol = chunkerkerneval(chnkr, ikern,dens_comb, targets);

targinfo = [];
targinfo.r = targets;
utrue = kern1(srcinfo, targinfo);


% plot 

figure(1)
clf

t = tiledlayout('flow');

nexttile
zztarg = NaN*zeros(size(xxtarg));
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
zztarg = NaN*zeros(size(xxtarg));
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
zztarg = NaN*zeros(size(xxtarg));
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


