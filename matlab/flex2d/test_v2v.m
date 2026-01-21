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

% figure(1); clf
% plot(S,rand(S.npatches,1))
% hold on
% plot(chnkr,'x-')
% view(0,90)

% V = eval_gauss(S.r);


% modified flexural (a \Delta^2 u - b \Delta u - c u) 

a = 1.1;
b = 0.7;
c = 1/pi;


nps = S.npts;
src = [S.r(1,3600); S.r(2,3600)];





zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));

zks = [zk1 zk2];


v2v = flex2d.v2v_matgen(S,zks,1e-8);
f = eval_gauss(S.r);
figure(1); clf
scatter(S.r(1,:),S.r(2,:),8,log10(f)); 
colorbar

val = v2v*f;
val = val(3600)



fun = @(r,th) adap_v2v(r,th,src,zks);
val_adap = integral2(fun,0,1,0,2*pi)





function val = eval_gauss(targ)

val = exp( - 20*targ(1,:).^2 - 20*targ(2,:).^2 );
val = val(:);

end

function val = adap_v2v(r,th,src,zks)
sz = size(r);    
r = r(:);
th = th(:);

xs = r.*cos(th);
ys = r.*sin(th);




targinfo = [];
targinfo.r = [xs; ys];


f = eval_gauss(targinfo.r);

srcinfo = [];
srcinfo.r = [src(1); src(2)];

val = chnk.flex2d.kern(zks,srcinfo,targinfo,'s');
val = val(:);
val = val.*r.*f;
val = reshape(val,sz);


end