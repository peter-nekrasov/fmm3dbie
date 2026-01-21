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


% zks = [1,2];

% val = chnk.flex2d.kern(zks,srcinfo,targinfo,'s');

% val

zk = 10;
val = chnk.flex2d.helmdiffgreen(zk,srcinfo.r,targinfo.r);
