src = zeros(2,1);
targ = zeros(2,1);


src(1) = 1.2;
src(2) = 2.3;

targ(1) = 2.5;
targ(2) = 3.5; 


nvec = [0.3;-0.2];
[val,grad] = chnk.flex2d.bhgreen(src, targ);
grad = squeeze(grad);


sum(grad.*nvec)