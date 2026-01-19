run('/Users/yuguan/software/chunkie/startup.m')

src = zeros(2,1);
targ = zeros(2,1);


src(1) = 1.2;
src(2) = 2.3;

targ(1) = 2.5;
targ(2) = 3.5; 


nvec = [0.3;-0.2];
[val,grad,hess] = chnk.flex2d.bhgreen(src, targ);
grad = squeeze(grad);

val 
sum(grad.*nvec)

nu = .3;
nu*(hess(:,1)+hess(:,3))+(1-nu)*(nvec(1,:).^2.'.*hess(:,1)+...
   2*(nvec(1,:).*nvec(2,:)).'.*hess(:,2)+nvec(2,:).^2.'.*hess(:,3))






