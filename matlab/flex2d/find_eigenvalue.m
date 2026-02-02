clear 
clc 






n = 50;
k = (0:n).';
x = 10*cos(pi*k/n)+10;


f = @(x) x.^2-.5;
vals = f(x);
vals = flipud(vals); 

p_cheb = chebfun(vals, [0, 20], 'chebkind', 2);
r_real = roots(p_cheb) 