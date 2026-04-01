S = geometries.disk([],[],[3 3 3],8);

nu = 0.3;

thetas = pi/24;

% thetas = 0;
targinfo = []; targinfo.r = [cos(thetas); sin(thetas); 0];
targinfo.n = targinfo.r;
targinfo.du = [-sin(thetas); cos(thetas); 0];
targinfo.kappa = 1;

Av2b = bh2d.v2b_matgen_neu(S,0,targinfo,1e-8);

h = 1e-6;
targinfo2 = []; targinfo2.r = targinfo.r + h*(0:2).*targinfo.n;

Aeval = bh2d.v2b_matgen_dir(S,0,targinfo2,1e-8);


dens = test_fn(S.r(1,:),S.r(2,:)).';

val = Av2b*dens;

d1 = [-3/2	2	-1/2]/h;
val2 = d1*(Aeval*dens);

err = abs(val - val2);

fprintf('neu relative error: %5.5e\n', vecnorm(err) / vecnorm(val));

bdry_pt = [cos(thetas); sin(thetas)];
bdry_n = [cos(thetas); sin(thetas)];
bdry_tau = [-sin(thetas); cos(thetas)];

h = 0.01;

d2dn2 = [	15/4	-77/6	107/6	-13	61/12	-5/6 0] / h^2;
d2dtau2 = [-1/12	4/3	-5/2	4/3	-1/12] / h^2;

d3dn3 = [ -49/8	29	-461/8	62	-307/8	13	-15/8 ]/h^3;
ddn = [-25/12	4	-3	4/3	-1/4 0 0]/h;

theta_t = atan2(bdry_n(2),bdry_n(1))-pi/2;
[xpts,ypts] = meshgrid(-2*h:h:2*h,0:h:6*h);
R = [cos(theta_t) -sin(theta_t); sin(theta_t) cos(theta_t)];
rot_pts = R*[xpts(:) ypts(:)].';
eval_pts = rot_pts + bdry_pt;

bctarginfo = []; bctarginfo.r = [eval_pts; eval_pts(1,:)*0];

Av2bsupp = bh2d.v2b_matgen_supp2(S,0,nu,targinfo,eps);

Aeval = bh2d.v2b_matgen_dir(S,0,bctarginfo,1e-8);

val = Av2bsupp*dens;

us = Aeval*dens;
us = reshape(us,size(xpts));

val2 = d2dn2*us(:,3) + nu*d2dtau2*us(1,:).';

err = abs(val -val2);

fprintf('supp2 relative error: %5.5e\n', vecnorm(err) / vecnorm(val));

Av2bfree = bh2d.v2b_matgen_free2(S,0,nu,targinfo,eps);

val = Av2bfree*dens;

d3dndtau2 = d2dtau2.*ddn.';

kappa =1;
val2 = d3dn3*us(:,3) + (2-nu)*sum(d3dndtau2.*us,'all') ...
   + (1-nu)*kappa*(d2dtau2*us(1,:).' - d2dn2*us(:,3));

err = abs(val-val2);

fprintf('free2 relative error: %5.5e\n', vecnorm(err) / vecnorm(val));

function val = test_fn(x,y)

    val = exp( - 100*x.^2 - 50*y.^4 );

end

