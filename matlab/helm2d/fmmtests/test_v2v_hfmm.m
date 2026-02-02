S = geometries.disk([],[],[4 4 4],8);

zk = pi;
%figure(1); clf
%plot(S,rand(S.npatches,1))

%%

test_fn = exp(-5*((S.r(1,:) - 0.3).^2+(S.r(2,:)+1/pi).^2));
test_fn = test_fn.';

tic;
A = helm2d.slp_matgen(S,zk,1e-12);
toc;

sol1 = A*test_fn;

tic;
[v2v_cor,nover] = helm2d.get_quad_cor_sub(S,zk, 1e-12);
toc;

sol2 = helm2d.apply_v2v(S,zk,test_fn,v2v_cor,nover,1e-12);

vecnorm(sol1 - sol2) / vecnorm(sol1)


