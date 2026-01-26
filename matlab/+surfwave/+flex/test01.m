%% Finite difference test for Green's functions
% Checking that the Green's functions satisfy the following relation:
%       1/2(\alpha \Delta^2 + \beta)G_S + \gamma G_\phi = 0
% away from the diagonal (x \neq y)

h = 0.05;
xs = 1:h:(1+8*h);
[X,Y] = meshgrid(xs);
alpha = 0.5;
beta = -1.5;
gamma = -2.5;
R = sqrt(X.^2 + Y.^2);

src = [];
src.r = [0;0];
targ = [];
targ.r = [X(:) Y(:)].';

[rts,ejs] = surfwave.flex.find_roots(alpha,beta,gamma);
[gs,~,~,~,fourth] = surfwave.flex.gsflex(rts,ejs,src,targ);
gphi = surfwave.flex.gphiflex(rts,ejs,src,targ);
bilapgs = fourth(:,:,1) + 2*fourth(:,:,3) + fourth(:,:,5);

gs = reshape(gs,size(X));
bilapgs = reshape(bilapgs,size(X));
gphi = reshape(gphi,size(X));

% finite difference stencils
d4 = zeros(9,1);
d4(1) = 7/240;	
d4(2) = -2/5;
d4(3) = 169/60;
d4(4) = -122/15;
d4(5) = 91/8;
d4(6) = -122/15;
d4(7) = 169/60;
d4(8) = -2/5;
d4(9) = 7/240;
d4 = d4/h^4; 

d2 = zeros(9, 1);
d2(1) = -1/560;
d2(2) = 8/315;
d2(3) = -1/5;
d2(4) = 8/5;
d2(5) = -205/72;
d2(6) = 8/5;
d2(7) = -1/5;
d2(8) = 8/315;
d2(9) = -1/560;
d2 = d2 / h^2;

bilap = zeros(9,9);
bilap(:,5) = d4;
bilap(5,:) = bilap(5,:) + d4.';
bilap = bilap + 2*(d2*d2.');

err = abs(0.5*alpha*sum(bilap.*gs,'all') - 0.5*beta*gs(5,5) + gamma*gphi(5,5)) / max(abs(gs(:))) 
err2 = abs(0.5*alpha*bilapgs(5,5) - 0.5*beta*gs(5,5) + gamma*gphi(5,5)) / max(abs(gs(:)))

%% test ifr2logr function

gsr2logr = surfwave.flex.gsflex(rts,ejs,src,targ,true);
gsr2logr = reshape(gsr2logr,size(X));

err3 = sum(bilap.*(gs - gsr2logr),'all')

xs = -4:0.01:4;
ys = 0*xs;

src = []; 
src.r = [0.1;0.2];
targ = [];
targ.r = [xs(:) ys(:)].';

rs = sqrt((src.r(1,:) - targ.r(1,:)).^2 + (src.r(2,:) - targ.r(2,:)).^2); 

gs = surfwave.flex.gsflex(rts,ejs,src,targ,false).';
gsr2logr = surfwave.flex.gsflex(rts,ejs,src,targ,true).'+1/(8*pi*alpha)*rs.^2.*log(abs(rs).^2);

err4 = max(gs - gsr2logr)


