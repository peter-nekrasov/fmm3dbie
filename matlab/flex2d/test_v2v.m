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



% modified flexural (a \Delta^2 u - b \Delta u - c u) 

a = 1.1;
b = 0.7;
c = 1/pi;


nps = S.npts;




zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));
zk1
zk2 
zks = [zk1 zk2];


v2v = flex2d.v2v_matgen(S,zks,1e-12);
f = eval_gauss(S.r);

% quadrature 
val_quad = v2v*f;


%% adaptive integration 
val_adap = zeros(size(val_quad));
for i=1:numel(val_adap)
  src = [S.r(1,i); S.r(2,i)];
  fun = @(r,th) adap_v2v(r,th,src,zks);
  val_adap(i) = integral2(fun,0,1,0,2*pi,"AbsTol",1e-12,"RelTol",1e-12);
end


rel_err = norm(val_quad-val_adap)/norm(val_adap)

%%


figure(2); clf
scatter(S.r(1,:),S.r(2,:),8,log10(abs(val_quad-val_adap)./abs(val_adap))); 
colorbar


%%


function val = eval_gauss(targ)
c1 = 0; 
c2 = 0;
val = exp( - 10*(targ(1,:)-c1).^2 - 10*(targ(2,:)-c2).^2 );
val = val(:);

end

function val = adap_v2v(r,th,src,zks)
% sz = size(r);    
% r = r(:);
% th = th(:);
% 
% xs = r.*cos(th);
% ys = r.*sin(th);
% 
% 
% 
% 
% targinfo = [];
% targinfo.r = [xs; ys];
% 
% 
% f = eval_gauss(targinfo.r);
% 
% srcinfo = [];
% srcinfo.r = [src(1); src(2)];
% 
% val = chnk.flex2d.kern(zks,srcinfo,targinfo,'s');
% val = val(:);
% val = val.*r.*f;
% val = reshape(val,sz);



xs = r .* cos(th);
ys = r .* sin(th);


sz = size(xs);


targinfo = [];
targinfo.r = [xs(:).'; ys(:).'];


f = eval_gauss(targinfo.r);   

srcinfo = [];
srcinfo.r = [src(1); src(2)];


k = chnk.flex2d.kern(zks,srcinfo,targinfo,'s');
k = k(:);


R = hypot(xs,ys);

val = k .* f .* R(:);

val = reshape(val,sz);


end