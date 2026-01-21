run('/Users/yuguan/software/chunkie/startup.m')

src = zeros(2,1);
targ = zeros(2,1);


src(1) = 1.2;
src(2) = 2.3;

targ(1) = 2;
targ(2) = 3.5; 


srcinfo = [];
srcinfo.r = src;


traginfo = [];
targinfo.r = targ;


zks = [1,0.002];

% val = chnk.flex2d.kern(zks,srcinfo,targinfo,'s');
% 
% val

zk1 = 10;
val1 = chnk.flex2d.helmdiffgreen(zk1,srcinfo.r,targinfo.r);


zk2 = 0.002;
val2 = chnk.flex2d.helmdiffgreen(zk2,srcinfo.r,targinfo.r);
val1 
val2

val = 1/(zk1^2-zk2^2)*(val1-val2)