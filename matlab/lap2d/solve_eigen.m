
clear 
clc 


run('/Users/yuguan/software/chunkie/startup.m')
run('/Users/yuguan/Dropbox/fmm3dbie/matlab/startup.m')
addpath '/Users/yuguan/Dropbox/fmm3dbie/src'


S = geometries.disk([],[],[4 4 4],8);
nch = 4*4;
cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = sort(chnkr);



g = sin(0.62*chnkr.r(1,:)) + ...
    cos(2.4*chnkr.r(1,:)) + ...
    cos(0.85*chnkr.r(1,:)).*sin(2.1*chnkr.r(2,:)) + ...
     cos(1.9*chnkr.r(1,:)).*cos(0.33*chnkr.r(2,:)) + ...
     sin(-2.6*chnkr.r(1,:)).*sin(1.45*chnkr.r(2,:))+1;
g = g(:);

f = cos(1.75*chnkr.r(2,:)) + ...
    sin(1.3*chnkr.r(1,:)).*cos(0.72*chnkr.r(2,:)) + ...
     cos(-2.2*chnkr.r(1,:)).*sin(0.58*chnkr.r(2,:)) + ...
     cos(3.1*chnkr.r(1,:)).*cos(-0.41*chnkr.r(2,:))+1;
f = f(:);



n = 20;
k = (0:n).';

eps0 = 1e-4;
a = 4;
b = 8;
x = (a+b)/2+(b-a)/2*cos(pi*(n-k)/n);

vals = zeros(n+1,1);

for i=1:n+1

    h2d_d = @(s,t) chnk.helm2d.kern(sqrt(x(i)), s, t, 'd');
    sysmat = -0.5*eye(chnkr.npt)+chunkermat(chnkr, h2d_d); 
    
    
    h = sysmat\g;
    vals(i) = 1/sum(h.*f.*chnkr.wts(:));

end


% p_cheb = chebfun(vals, [a,b], 'chebkind', 2);
% tol = 1e-7;
% q = real(p_cheb.*conj(p_cheb));
% dq = diff(q);
% xc = roots(dq);
% xc = [a; xc; b];
% xc = xc(xc>=a & xc<=b);
% qc = q(xc);
% idx = qc<tol^2;
% rs = xc(idx);
% q(rs)
% rs
% sqrt(rs)


%%

figure 
plot(x,real(vals))
hold on 
plot(x,imag(vals))
legend('real','imag')

%%

aa = -1;
bb = 1;
x = (aa+bb)/2+(bb-aa)/2*cos(pi*(n-k)/n);
[N,X] = ndgrid(k,x);
c2v = cos(N.*acos(X)).';
% v2c = inv(c2v);
cfs = c2v\vals;

%%

coll = zeros(n,n);
for ii=1:(n-1)
    coll(ii,ii+1) = 1/2;
    coll(ii+1,ii) = 1/2;
end

coll(1,2) = 1;

cfred = cfs(1:(end-1))/(2*cfs(end));

coll(end,:)=coll(end,:)-cfred.';

eis = eigs(coll,n);

eis = eis(abs(real(eis))<1);
tol = 1E-7;
eis = eis(abs(imag(eis))<tol);


(eis-aa)/2*(b-a)+a

%%

% zk = 2.404825557;
% h2d_d = @(s,t) chnk.helm2d.kern(zk, s, t, 'd');
% sysmat = -0.5*eye(chnkr.npt)+chunkermat(chnkr, h2d_d); 
% cond(sysmat)
% h = sysmat\g;
% 1/sum(h.*f)
