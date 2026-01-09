%% Get chnkr objects

verts_all = load('geo.mat');
verts = verts_all.pts.';
verts = (verts);
% verts = verts(:,1:12:end-100);
verts = verts(:,1:30:end-160);
verts(2,:) = 0.5*verts(2,:);
nchs = ones(size(verts,2),1); 

opts = [];
opts.nchs = 2;
opts.lam  = 3;
opts.step_fact = 0.9;
opts.n_newton = 500;
opts.etol = 1E-14;
[chnkr,err] = chnk.smoother.smooth(verts,opts);

chnkr.npt 

figure(1); clf 
t = tiledlayout(1,3);
nexttile
plot(chnkr,'k-')
hold on
quiver(chnkr)

opts = []; opts.mv_bdries = 1;
[chnkr1,eps] = arcresample(chnkr,opts);
chnkr1 = refine(chnkr1,struct('maxchunklen',0.02));
chnkr1 = sort(chnkr1);
chnkr1.npt 

opts = []; opts.mv_bdries = 1;
[chnkr2,eps] = arcresampler(chnkr,opts);
chnkr2 = refine(chnkr2,struct('maxchunklen',0.02));
chnkr2 = sort(chnkr2);

% kappa = signed_curvature(chnkr);
% kappamax = max(abs(kappa)).*(chunklen(chnkr)).';
% [~,idx] = maxk(kappamax,floor(0.4*chnkr.nch));
% figure(1); clf 
% hold on
% plot(chnkr,'k-')
% hold on
% quiver(chnkr)
% scatter(chnkr.r) 


% opts = [];
% % opts.maxchunklen = 10;
% % opts.splitchunks = idx;
% % opts.lvlr = 'n';
% chnkr = refine(chnkr, opts);
% 
% chnkr.npt


% verts2 = verts_all.verts2;
% nchs = ones(size(verts2,2),1); 
% 
% opts = [];
% opts.nchs = nchs;
% opts.lam  = 5;
% opts.step_fact = 0.9;
% opts.n_newton = 500;
% opts.etol = 1E-6;
% [chnkr2,err] = chnk.smoother.smooth(verts2,opts);


nexttile
plot(chnkr1,'k-')
hold on
quiver(chnkr1)
nexttile
plot(chnkr2,'k-')
hold on
quiver(chnkr2)
% 
% chnkr = chnkr1;
% 
% figure(2); clf 
% [xleg,wleg,v2c,~] = lege.exps(16);
% coefs = v2c*squeeze(chnkr.r(1,:,:));
% plot(log10(abs(coefs)))

%%

k = 16;

rt = chnkr1.r;
plot(rt(1,:), rt(2,:), 'k.')

x = rt(1,:);
y = rt(2,:);


x = reshape(x,k,[]);
y = reshape(y,k,[]);


[xl,wl,ul,vl] = lege.exps(k);
ps = lege.pols(-1,k-1);
vx = ps.'*ul*x;
vy = ps.'*ul*y;

[rs,~] = chunkends(chnkr1);
rv = squeeze(rs(:,1,:));
%%

pgon = polyshape(vx, ...
    vy);
tr = triangulation(pgon);
gm = fegeometry(tr);
gm = addVertex(gm,"Coordinates",rv.');
gm = generateMesh(gm,'GeometricOrder','linear','Hmax',0.03,'Hgrad',1.2);

pdemesh(gm); hold on;
msh = gm.Mesh;

% find all boundary edges

ee = msh.Elements;
nn = msh.Nodes;
tr_tmp = triangulation(ee.', nn.');

ff = freeBoundary(tr_tmp).';

% scatter(msh.Nodes(1,ff(1,1:(end-1))),msh.Nodes(2,ff(1,1:(end-1))),...
%     20,'filled');
% 
% scatter(vx,vy,20,'*');

node_vert = findNodes(msh,"nearest",rv);

nv = numel(node_vert);

%ff = fliplr(flipud(ff));
n1 = node_vert(1);

ind1 = find(ff(1,:)==n1);
ff = circshift(ff,[0,-ind1+1]);

[I,J] = meshgrid(node_vert,ff(1,:));
[inode,iff] = find(I.'==J.');
[inode,iord ] = sort(inode);
iff = iff(iord);

in_st_ou = node_vert(1);
in_fi_ou = node_vert(2);
%%
node_vert = [node_vert,node_vert(1)];
icount = 1;

xs = zeros(k,size(ff,2));
ys = zeros(k,size(ff,2));
iverts = zeros(2,size(ff,2));
nbedge = size(ff,2);



strts = zeros(nbedge,1);
finis = strts;
iis   = strts;

for ii=1:(nv)
    i0 = iff(ii);
    if (ii ~= nv)
        i1 = iff(ii+1);
        inds = ff(1,i0:(i1-1));
    else
        i1 = iff(1);
        inds = ff(1,i0:end);
    end
    v1 = nn(:,ff(1,i0));
    v2 = nn(:,ff(1,i1));
    vdist = norm(v2-v1);
    istrts = nn(:,inds);
    iends  = nn(:,[inds(2:end),ff(1,i1)]);
    sdist1 = vecnorm(istrts-v1);
    sdist2 = vecnorm(istrts-v2);
    spar = sdist1./(sdist1+sdist2);

    fdist1 = vecnorm(iends-v1);
    fdist2 = vecnorm(iends-v2);
    fpar = fdist1./(fdist1+fdist2);

    for iii=2:numel(inds)
        xloc = x(:,ii);
        yloc = y(:,ii);
        tt = spar(iii)*2-1;
        xm = (lege.pols(tt,k-1)).'*ul*xloc;
        ym = (lege.pols(tt,k-1)).'*ul*yloc;
        nn(:,inds(iii)) = [xm;ym];
    end

    inds = [inds,ff(1,i1)];
    for iii=1:numel(inds)-1
        xloc = x(:,ii);
        yloc = y(:,ii);
        strt = 2*spar(iii)-1;
        strts(icount) = strt;
        fini = 2*fpar(iii)-1;
        finis(icount) = fini;
        iis(icount)   = ii;
        ts = (xl+1)/2*(fini-strt)+strt;
        xts = (lege.pols(ts,k-1)).'*ul*xloc;
        yts = (lege.pols(ts,k-1)).'*ul*yloc;
        xs(:,icount) = xts;
        ys(:,icount) = yts;
        iverts(1,icount) = inds(iii);
        iverts(2,icount) = inds(iii+1);
        icount = icount + 1;
    end
end


%%

norder = 20;
[uvs] = koorn.rv_nodes(norder);
[wts] = koorn.rv_weights(norder);
nuvs = size(uvs,2);

nelem = size(ee,2);

xx = zeros(nuvs,nelem);
yy = zeros(nuvs,nelem);

edgepatchinfo = [];
point_id = [];
patch_id = [];
u_vals = [];

for ii=1:nelem

    nds = ee(:,ii);
    n1 = nds(1);
    n2 = nds(2);
    n3 = nds(3);
    v1 = nn(:,n1);
    v2 = nn(:,n2);
    v3 = nn(:,n3);

    [irw,icl] = find((iverts == (n1)) + (iverts == (n2)) + (iverts == (n3)));
    [i1,i2] = find((icl==icl.')-eye(numel(icl)));

    if (numel(i1)>0)

        enum = icl(i1(1));
        inot = setdiff([n1,n2,n3],iverts(:,enum));
        v1 = nn(:,inot);
        v2 = nn(:,iverts(1,enum));
        v3 = nn(:,iverts(2,enum));

        xt = xs(:,enum);
        yt = ys(:,enum);

        strt = strts(enum);
        fini = finis(enum);
        iiv  = iis(enum);

        edgepatchinfo = [edgepatchinfo,[strt,fini,iiv,ii].'];
        xl = lege.exps(16);
        inds = find((xl >= strt).*(xl<= fini));
        is_in_patch = find((xl >= strt).*(xl<= fini));
        xloc = xl(is_in_patch);
        %xrel = ((xloc-strt)/(fini-strt)*2-1);
        xrel = 1-(xloc-strt)/(fini-strt);
        u_vals  = [u_vals;xrel];
        patch_id = [patch_id;ii*ones(size(xrel))];
        point_id = [point_id,inds.' + (iiv-1)*16];

        % gamma(0)

        gam0 = [(lege.pols(-1,k-1)).'*ul*xt;(lege.pols(-1,k-1)).'*ul*yt];
        gamL = [(lege.pols( 1,k-1)).'*ul*xt;(lege.pols( 1,k-1)).'*ul*yt];
        gamx = [(lege.pols((1-uvs(1,:).')*2-1,k-1)).'*ul*xt,...
                (lege.pols((1-uvs(1,:).')*2-1,k-1)).'*ul*yt].';

        ntrm = 1 - uvs(1,:)-uvs(2,:);
        xi   = uvs(1,:);
        eta  = uvs(2,:);
        dtrm = 1 - uvs(1,:);
        tran = ntrm.*gamL+xi.*gam0 + eta.*v1 + ...
            ntrm./dtrm.*(gamx-dtrm.*gamL-xi.*gam0);
        xx(:,ii) = tran(1,:);
        yy(:,ii) = tran(2,:);


    else

        uvec = v3-v1;
        vvec = v2-v1;
        xvec = v1(1) + uvec(1)*uvs(1,:) + vvec(1)*uvs(2,:);
        yvec = v1(2) + uvec(2)*uvs(1,:) + vvec(2)*uvs(2,:);   
        xx(:,ii) = xvec;
        yy(:,ii) = yvec;

    end

end

v_vals = 0*u_vals;
uv_bndry = [u_vals.';v_vals.'];

amat = koorn.vals2coefs(norder,uvs);
cfx  = amat*xx;
cfy  = amat*yy;

nord = 8;
[uv] = koorn.rv_nodes(nord);
[wt] = koorn.rv_weights(nord);
[pols,dersu,dersv] = koorn.ders(norder, uv);

srcx   = pols.'*cfx;
srcy   = pols.'*cfy;
srcdxu = dersu.'*cfx;
srcdxv = dersv.'*cfx;
srcdyu = dersu.'*cfy;
srcdyv = dersv.'*cfy;
srcdet = diag(wt)*(srcdxu.*srcdyv - srcdxv.*srcdyu);

npts = numel(srcdet);

srcvals = zeros(12,npts);
srcvals(1,:) = srcx(:);
srcvals(2,:) = srcy(:);
srcvals(4,:) = srcdxu(:);
srcvals(5,:) = srcdyu(:);
srcvals(7,:) = srcdxv(:);
srcvals(8,:) = srcdyv(:);
srcvals(12,:) = -1;


srfr = surfer(nelem,nord,srcvals,1);
errs = surf_fun_error(srfr,srcvals(4,:));

%%

figure(4); 
plot(srfr,rand(srfr.npatches,1))
hold on
plot(chnkr1,'x-')
hold on
quiver(chnkr1)
view(0,90)

chnkr = chnkr1;
S = srfr;

[~,I] = sort(point_id);

uv_bndry = uv_bndry(:,I);
patch_id = patch_id(I);

return

save("spiky_blob_geom.mat","chnkr","S","uv_bndry","patch_id")