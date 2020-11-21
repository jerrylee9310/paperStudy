function m = invgridkb(d,k,n,osf,wg,beta)

% function m = gridkb(d,k,w,n,osf,kw,opt)
%
%     d -- cartesian k-space data
%     k -- k-trajectory, scaled -0.5 to 0.5
%     n -- image size (m will be osf*n X osf*n)
%     osf -- oversampling factor (usually between 1 and 2)
%     wg -- full kernel width in oversampled grid samples (usually 3 to 7)
%
%     m -- gridded k-space data (radial)
%
%  Uses optimum Kaiser-Bessel window for a given
%    oversampling factor and kernel size
%  Now uses Phil's numbers
[fov,spk] = size(k);

% convert to single column
k = k(:);

% width of the kernel on the original grid
% kw = wg*osf;
kw = wg/osf;

% compute kernel, assume e1 is 0.001, assuming nearest neighbor
kosf = floor(0.91/(osf*1e-3));

% half width in undersampled grid units
% kwidth = kw/(osf*2);
kwidth = osf*kw/2;

% beta from the Beatty paper
% beta = pi*sqrt((kw*(osf-0.5)).^2-0.8);

% compute kernel
om = [0:kosf*kwidth]/(kosf*kwidth);
p = besseli(0,beta*sqrt(1-om.*om));
p = p./p(1);
% last sample is zero so we can use min() below for samples bigger than kwidth
p(end) = 0;

% convert k-space samples to matrix indices
nx = (n*osf/2+1) + osf*n*real(k);
ny = (n*osf/2+1) + osf*n*imag(k);

rx = [1:fov].';
ry = [1:spk];

radn = rx + 1j*ry;
radx = real(radn);
rady = imag(radn);
radx = radx(:);
rady = rady(:);

m = zeros(fov,spk);

% loop over samples in kernel at grid spacing
for lx = -kwidth:kwidth,
  for ly = -kwidth:kwidth,

    % find nearest samples
    nxt = round(nx+lx);
    nyt = round(ny+ly);

    % seperable kernel value
    kkx = min(round(kosf*abs(nx-nxt)+1), floor(kosf*kwidth)+1);
    kwx = p(kkx);
    kky = min(round(kosf*abs(ny-nyt)+1), floor(kosf*kwidth)+1);
    kwy = p(kky);

    % if data falls outside matrix, put it at the edge, zero out below
    nxt = max(nxt,1); nxt = min(nxt,osf*n);
    nyt = max(nyt,1); nyt = min(nyt,osf*n);

    % accumulate gridded data
%     m = m+sparse(nxt,nyt,dw.*kwx'.*kwy',osf*n,osf*n);
    for nd = 1:length(k)
        m(radx(nd),rady(nd)) = m(radx(nd),rady(nd)) + d(nxt(nd),nyt(nd))*kwx(nd)*kwy(nd);
    end
    
  end
end

