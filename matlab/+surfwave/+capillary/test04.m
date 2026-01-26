% Check that S[G_S] = G_phi

alpha = 2; beta = -1+0.8i; gamma = 2;
[rts,ejs] = surfwave.capillary.find_roots(alpha,beta,gamma);

% G_S centered at [1;1]

x = [3;2]; % target
xp = [1;1]; % source

gp = surfwave.capillary.gphihelm(rts,ejs,xp,x);

integrand = @(r,t) r.*gs(rts,ejs,xp,r,t)./sqrt( (r.*cos(t) - x(1)).^2 + (r.*sin(t) - x(2)).^2) / (4*pi);

rmin = 0;
rmax = 1000;
tmin = 0;
tmax = 2*pi;

result = integral2(integrand, rmin, rmax, tmin, tmax,"AbsTol",1e-8);

err = abs(result - gp) / abs(gp)

function val = gs(rts,ejs,xp,r,t)

    targs = [r(:).*cos(t(:)), r(:).*sin(t(:))].';
    gs = surfwave.capillary.gshelm(rts,ejs,xp,targs);
    val = reshape(gs, size(r));
    
end