%% Get chnkr objects

verts_all = load('geo.mat');
verts = verts_all.pts.';
verts = (verts);
% verts = verts(:,1:12:end-100);
verts = verts(:,1:30:end-160);
verts(2,:) = 0.5*verts(2,:);
nchs = ones(size(verts,2),1); 

opts = [];
opts.nchs = 2;
opts.lam  = 3;
opts.step_fact = 0.9;
opts.n_newton = 500;
opts.etol = 1E-14;
[chnkr,err] = chnk.smoother.smooth(verts,opts);

chnkr.npt 

figure(1); clf 
plot(chnkr,'k-')
hold on
quiver(chnkr)
hold on 
plot(verts(1,:),verts(2,:),'kx')
