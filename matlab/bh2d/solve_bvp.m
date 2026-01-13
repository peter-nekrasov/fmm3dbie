% genpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/chunkie/')

clear 
clc

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

figure(1); clf
plot(S,rand(S.npatches,1))
hold on
plot(chnkr,'x-')
view(0,90)

V = eval_gauss(S.r);

zk = 0;

%% v2v
% A = lap2d.slp_matgen(S,1e-9);
% lhs_11 = -eye(S.npts) + V.*A;

A = bh2d.v2v_matgen(S,zk,1e-8);
lhs_11 = -eye(S.npts) + V.*A;


% val = A*V;
% i = 3600;
% 
% xt = S.r(1,i);
% yt = S.r(2,i);
% 
% 
% val_quad = val(i);
% val_adap = v2v_adap(xt,yt);
% abs(val_quad-val_adap)/abs(val_adap)





%% b2v


fkern = @(s,t) bh2d.kern(zk,s,t,'s');
lhs_12 = V.*chunkerkernevalmat(chnkr,fkern,S.r(1:2,:));



% fkern_lap = @(s,t) chnk.lap2d.kern(s,t,'s');
% A = chunkerkernevalmat(chnkr,fkern,S.r(1:2,:));
% src = chnkr.r(1:2,:);





%% v2b
targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
lhs_21 = bh2d.gdn_matgen(S,zk,targinfo,1e-8);



%% b2b
% l2d_sp = kernel('l','sp');
% lhs_22 = 0.5*eye(chnkr.npt)+chunkermat(chnkr,l2d_sp); %0.5*eye(n)... -

bh2d_sp = @(s,t) bh2d.kern(zk,s,t,'sp');
lhs_22 = 0.5*eye(chnkr.npt)+chunkermat(chnkr,bh2d_sp);





%%

lhs = [lhs_11, lhs_12; lhs_21, lhs_22];
% 
% 
% 
rhs_1 = (4+V(:).').*sin(S.r(1,:)).*sin(S.r(2,:)) ; 
rhs_2 = cos(chnkr.r(1,:)).*sin(chnkr.r(2,:)).*chnkr.n(1,:) ...
          + sin(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(2,:) ; % sin(xs(1,:))

rhs = [rhs_1 rhs_2].';

dens = lhs\rhs;

mu = dens(1:S.npts);
rho = dens(S.npts+1:end);

fkern_s = kernel('l','s');
u = A*mu + chunkerkerneval(chnkr,fkern_s,rho,S.r(1:2,:));

ref_u = sin(S.r(1,:)).*sin(S.r(2,:));
err = abs(u - ref_u(:)) / max(abs(u));

figure(1); clf
scatter(S.r(1,:),S.r(2,:),8,log10(err)); colorbar

fprintf('relative L^2/L^2 error: %5.5e\n', vecnorm(err .* S.wts(:)) / vecnorm(u .* S.wts(:)) );

return



%%
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end