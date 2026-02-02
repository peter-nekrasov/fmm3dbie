S = geometries.disk([],[],[4 4 4],8);

zk = pi;
% chnkr = chunkerfunc(@(t) starfish(t,5,0,[0,0],0,1));
nch = 4*4;
cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = sort(chnkr);

h2d_s = kernel('h','s',zk);

Ab2v = chunkerkernevalmat(chnkr,h2d_s,S.r(1:2,:));

test_fn = chnkr.r(1,:).*chnkr.r(2,:);
test_fn = test_fn(:);

sol1 = Ab2v * test_fn;

% A_native = l2d_s.eval(chnkr,S).*chnkr.wts(:).'; 
% sol2 = A_native*test_fn(:) + Ab2v_cor*test_fn;

opts.corrections = true;
Ab2v_cor = chunkerkernevalmat(chnkr,h2d_s,S.r(1:2,:),opts);

sol2 = helm2d.apply_b2v_neu(S,zk,chnkr,test_fn,Ab2v_cor,1e-12);

vecnorm(sol1 - sol2) / vecnorm(sol1)


h2d_d = kernel('h','d',zk);

Ab2v = chunkerkernevalmat(chnkr,h2d_d,S.r(1:2,:));

test_fn = chnkr.r(1,:).*chnkr.r(2,:);
test_fn = test_fn(:);

sol1 = Ab2v * test_fn;

% A_native = l2d_s.eval(chnkr,S).*chnkr.wts(:).'; 
% sol2 = A_native*test_fn(:) + Ab2v_cor*test_fn;

opts.corrections = true;
Ab2v_cor = chunkerkernevalmat(chnkr,h2d_d,S.r(1:2,:),opts);

sol2 = helm2d.apply_b2v_dir(S,zk,chnkr,test_fn,Ab2v_cor,1e-12);

vecnorm(sol1 - sol2) / vecnorm(sol1)
