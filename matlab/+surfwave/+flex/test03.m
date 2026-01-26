%% test for derivatives of greens functions
% part 1 - check helmholtz kernels
% 

k = 3;

src = [];
src.r = [5;0];

h = 0.01;

[X,Y] = meshgrid(0:h:8*h); 
targ = [];
targ.r = [X(:) Y(:)].';

d1 = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280]/h;
d2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560]/h^2;

[val,grad,hess] = surfwave.flex.helmdiffgreen(k,src.r,targ.r);
val = reshape(val,size(X));
gradx = reshape(grad(:,:,1),size(X));
grady = reshape(grad(:,:,2),size(X));
hessxx = reshape(hess(:,:,1),size(X));
hessxy = reshape(hess(:,:,2),size(X));
hessyy = reshape(hess(:,:,3),size(X));

err1 = d1*val(5,:).' - gradx(5,5);
err2 = d2*val(5,:).' - hessxx(5,5);

fprintf('error 1: %5.2e\n', abs(err1));
fprintf('error 2: %5.2e\n', abs(err2));


%% part 2 - nonlocal helmholtz kernels (G_S)

alpha = 2;
beta = 2;
gamma = -3;
[rts,ejs] = surfwave.flex.find_roots(alpha,beta,gamma);

src = [];
src.r = [1.1;0];

h = 0.01;

[X,Y] = meshgrid(0:h:8*h); 
targ = [];
targ.r = [X(:) Y(:)].';

[val] = surfwave.flex.gsflex(rts,ejs,src,targ);
[~,grad] = surfwave.flex.gsflex(rts,ejs,src,targ);
[~,~,hess] = surfwave.flex.gsflex(rts,ejs,src,targ);
[~,~,~,third] = surfwave.flex.gsflex(rts,ejs,src,targ);
[~,~,~,~,fourth] = surfwave.flex.gsflex(rts,ejs,src,targ);

val = reshape(val,size(X));
gradx = reshape(grad(:,:,1),size(X));
grady = reshape(grad(:,:,2),size(X));
hessxx = reshape(hess(:,:,1),size(X));
hessxy = reshape(hess(:,:,2),size(X));
hessyy = reshape(hess(:,:,3),size(X));
thirdxxx = reshape(third(:,:,1),size(X));
thirdxxy = reshape(third(:,:,2),size(X));
thirdxyy = reshape(third(:,:,3),size(X));
thirdyyy = reshape(third(:,:,4),size(X));
fourthxxxx = reshape(fourth(:,:,1),size(X));
fourthxxxy = reshape(fourth(:,:,2),size(X));
fourthxxyy = reshape(fourth(:,:,3),size(X));
fourthxyyy = reshape(fourth(:,:,4),size(X));
fourthyyyy = reshape(fourth(:,:,5),size(X));

d1 = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280]/h;
d2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560]/h^2;

err1 = d1*val(5,:).' - gradx(5,5);
err2 = d2*val(5,:).' - hessxx(5,5);

err3 = d1*val(:,5) - grady(5,5);
err4 = d2*val(:,5) - hessyy(5,5);

err5 = sum((d1.'*d1).*val,'all') - hessxy(5,5);

err6 = d1*hessxx(5,:).' - thirdxxx(5,5);
err7 = d1*hessxx(:,5) - thirdxxy(5,5);
err8 = d1*hessyy(5,:).' - thirdxyy(5,5);
err9 = d1*hessyy(:,5) - thirdyyy(5,5);

err10 = d1*thirdxxx(5,:).' - fourthxxxx(5,5);
err11 = d1*thirdxxy(5,:).' - fourthxxxy(5,5);
err12 = d1*thirdxyy(5,:).' - fourthxxyy(5,5);
err13 = d1*thirdyyy(5,:).' - fourthxyyy(5,5);
err14 = d1*thirdyyy(:,5) - fourthyyyy(5,5);

fprintf('error 1: %5.2e\n', abs(err1));
fprintf('error 2: %5.2e\n', abs(err2));
fprintf('error 3: %5.2e\n', abs(err3));
fprintf('error 4: %5.2e\n', abs(err4));
fprintf('error 5: %5.2e\n', abs(err5));
fprintf('error 6: %5.2e\n', abs(err6));
fprintf('error 7: %5.2e\n', abs(err7));
fprintf('error 8: %5.2e\n', abs(err8));
fprintf('error 9: %5.2e\n', abs(err9));
fprintf('error 10: %5.2e\n', abs(err10));
fprintf('error 11: %5.2e\n', abs(err11));
fprintf('error 12: %5.2e\n', abs(err12));
fprintf('error 13: %5.2e\n', abs(err13));
fprintf('error 14: %5.2e\n', abs(err14));



%% part 3 - nonlocal helmholtz kernels (G_\phi) 

alpha = 2;
beta = 2;
gamma = -3;
[rts,ejs] = surfwave.flex.find_roots(alpha,beta,gamma);

src = [];
src.r = [1.1;0];

h = 0.01;

[X,Y] = meshgrid(0:h:8*h); 
targ = [];
targ.r = [X(:) Y(:)].';

[val] = surfwave.flex.gphiflex(rts,ejs,src,targ);
[~,grad] = surfwave.flex.gphiflex(rts,ejs,src,targ);
[~,~,hess] = surfwave.flex.gphiflex(rts,ejs,src,targ);
[~,~,~,third] = surfwave.flex.gphiflex(rts,ejs,src,targ);
[~,~,~,~,bilap] = surfwave.flex.gphiflex(rts,ejs,src,targ);

val = reshape(val,size(X));
gradx = reshape(grad(:,:,1),size(X));
grady = reshape(grad(:,:,2),size(X));
hessxx = reshape(hess(:,:,1),size(X));
hessxy = reshape(hess(:,:,2),size(X));
hessyy = reshape(hess(:,:,3),size(X));
thirdxxx = reshape(third(:,:,1),size(X));
thirdxxy = reshape(third(:,:,2),size(X));
thirdxyy = reshape(third(:,:,3),size(X));
thirdyyy = reshape(third(:,:,4),size(X));
bilap = reshape(bilap,size(X));

d1 = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280]/h;
d2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560]/h^2;

err1 = d1*val(5,:).' - gradx(5,5);
err2 = d2*val(5,:).' - hessxx(5,5);

err3 = d1*val(:,5) - grady(5,5);
err4 = d2*val(:,5) - hessyy(5,5);

err5 = sum((d1.'*d1).*val,'all') - hessxy(5,5);

err6 = d1*hessxx(5,:).' - thirdxxx(5,5);

err7 = d1*hessxx(:,5) - thirdxxy(5,5);
err8 = d1*hessyy(5,:).' - thirdxyy(5,5);
err9 = d1*hessyy(:,5) - thirdyyy(5,5);

err10 = d2*hessxx(5,:).' + 2*d2*hessxx(:,5) + d2*hessyy(:,5) - bilap(5,5);


fprintf('error 1: %5.2e\n', abs(err1));
fprintf('error 2: %5.2e\n', abs(err2));
fprintf('error 3: %5.2e\n', abs(err3));
fprintf('error 4: %5.2e\n', abs(err4));
fprintf('error 5: %5.2e\n', abs(err5));
fprintf('error 6: %5.2e\n', abs(err6));
fprintf('error 7: %5.2e\n', abs(err7));
fprintf('error 8: %5.2e\n', abs(err8));
fprintf('error 9: %5.2e\n', abs(err9));
fprintf('error 10: %5.2e\n', abs(err10));