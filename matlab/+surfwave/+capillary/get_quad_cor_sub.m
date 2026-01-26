function xmat = get_quad_cor_sub(S, gs_kern, eps, zpars, type, norderup)

    if strcmp(type,'gs')
        iker = 0;
    elseif strcmp(type,'gphi')
        iker = 1;
    elseif strcmp(type,'lapgs')
        iker = 2;
    elseif strcmp(type,'lapgphi')
        iker = 3;
    elseif strcmp(type,'s3dgphi')
        iker = 5;
    else
        error('kernel name not recognized') 
    end 

    if nargin < 5
        ivpp = 1;
    end

    gs_kern = kernel(gs_kern);
    if nargin < 6
        norderup = 0;
    end
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    npatches = S.npatches;
    npts = S.npts;
    ndtarg = 3;

    % [patch_id, uvs_targ] = get_patch_id_uvs(S);

    %this might need fixing

    iptype_avg = floor(sum(iptype)/(npatches+0.0d0));
    norder_avg = floor(sum(norders)/(npatches+0.0d0));

    [rfac, rfac0] = get_rfacs(norder_avg,iptype_avg);
    % prin2('rfac = *',rfac,1)
    % prin2('rfac0 = *',rfac0,1)

        % allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

    [cms, rads] = get_centroid_rads(npatches,norders,ixyzs,iptype,npts, ... 
       srccoefs);

    rad_near = rads*rfac;

    % 
    % find near quadrature correction interactions
    % 

    targs = srcvals(1:3,:);
    ntarg = npts;

    % prinf('entering find near mem',0,0)
    nnz = findnearmem(cms, npatches, rad_near, ndtarg, targs, npts);
    % prinf('nnz = *',nnz,1);

    [row_ptr,col_ind] = findnear(cms,npatches,rad_near,ndtarg,targs,npts,nnz);

    iquad = get_iquad_rsc(npatches,ixyzs,npts,ntarg,nnz,row_ptr,col_ind);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wts = get_qwts_f(npatches,norders,ixyzs,iptype,ntarg,srcvals);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   compute near quadrature correction

    nquad = iquad(nnz+1)-1;

    iquadtype = 1;

    wnear = cell(length(iker),1);
    for i = 1:length(iker)
        if iker < 6
        wnear{i} = surfwave.capillary.getnearquad_capillary(npatches,norders,ixyzs, ...
              iptype,npts,srccoefs,srcvals,eps,iquadtype,nnz, ...
              row_ptr,col_ind,iquad,rfac0,zpars,nquad,iker(i)) ;
        elseif iker == 6
            opts.rep = 'eval';
            opts.format = 'rsc';
            wstruc = lap3d.neumann.get_quadrature_correction(S,eps,S,opts);
            wnear{i} = wstruc.wnear;
        end
    end

    if length(iker) == 1
        wnear = wnear{1};
    end

    npols = ixyzs(2:end)-ixyzs(1:end-1);
    npts_col = npols(col_ind);

    [nt,~] = size(row_ptr);
    ntarg = nt-1;
    nrep = zeros(ntarg,1);
    istarts = row_ptr(1:end-1);
    iends = row_ptr(2:end)-1;
    icol_ind = zeros(sum(npts_col),1);
    isrcinds = cell(npatches,1);
    for i=1:npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
    end

    istart = 1;
    for i=1:ntarg
        nrep(i) = sum(npts_col(istarts(i):iends(i))); 
        iinds = horzcat(isrcinds{col_ind(istarts(i):iends(i))});
        nelem = length(iinds);
        icol_ind(istart:istart+nelem-1) = iinds;
        istart = istart+nelem;
    end
    irow_ind = repelem((1:ntarg)',nrep);    
    
    xmat = sparse(irow_ind,icol_ind, wnear, ntarg, S.npts);

    Asmth_over = smooth_sparse_quad(gs_kern,targs,S,row_ptr,col_ind,norderup); 

    xmat = xmat - Asmth_over;

end

% tic;
%    Gsquad = conv_rsc_to_spmat(S,row_ptr,col_ind,wnear);
% 
%     [i,j,vals] = find(Gsquad);
% 
%     nnz = length(i);
%     lbat = 1e3;
%     nbat = ceil(nnz/lbat);
%     % tic;
% 
%     vals_sub = zeros(size(vals),'like',vals);
%     for k = 1:nbat
%         ks = (lbat*(k-1)+1):min(lbat*k,nnz);
%         rsrc = S.r(:,i(ks));
%         rtarg = S.r(:,j(ks));
%         rnear = rtarg - rsrc;
%         vals_sub(ks) = gs_kern(struct('r',[0;0;0]),struct('r',rnear)).*S.wts(j(ks));
%     end
%     % rsrc = S.r(:,i);
%     % rtarg = S.r(:,j);
%     % rnear = rtarg - rsrc;
%     % vals_sub = gs_kern(struct('r',[0;0;0]),struct('r',rnear)).*S.wts(j);
%     xmat2 = sparse(i,j, vals-vals_sub);
%     wnear2 = wnear;
% toc;
%
%
%
%
