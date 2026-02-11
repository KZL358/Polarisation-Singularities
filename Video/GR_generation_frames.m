%% ================================
%  PART 1 — PRECOMPUTE ALL FRAMES
%  ================================
clear; clc; close all;
%% Parameters
nWaves    = 40;
lambda    = 0.5;
kMag      = 2*pi/lambda;
gridN     = 100;
domainMin = -0.5;
domainMax =  0.5;
rng(2)
% Animation α values
nFrames   = 300;
ampList   = linspace(0,1,nFrames);
%% Build grid
x = linspace(domainMin, domainMax, gridN);
y = linspace(domainMin, domainMax, gridN);
z = linspace(domainMin, domainMax, gridN);
[X,Y,Z] = ndgrid(x,y,z);
%% Prepare GW fields
h_base = complex(zeros(gridN,gridN,gridN,3,3));
h_nth  = complex(zeros(gridN,gridN,gridN,3,3));
%% Generate random GW waves
fprintf("Generating random waves...\n");
A = [1 0 0]; % helper vector
for iWave = 1:nWaves
    k1 = [0 0 0];
    while norm(k1)==0
        k1 = randi([-10,10],[1 3]);
    end
    kVec = k1/norm(k1) * kMag;
    % Build TT basis
    if norm(cross(kVec,A))<1e-12
        Ause=[0 1 0];
    else
        Ause=A;
    end
    e1 = cross(kVec,Ause); e1=e1/norm(e1);
    e2 = cross(kVec,e1);   e2=e2/norm(e2);
    Ep = e1.'*e1 - e2.'*e2;
    Ex = e2.'*e1 + e1.'*e2;
    c1 = rand + 1i*rand;
    c2 = rand + 1i*rand;
    H  = c1*Ep + c2*Ex;
    phi = 2*pi*rand;
    phase    = kVec(1)*X + kVec(2)*Y + kVec(3)*Z + phi;
    expPhase = exp(1i*phase);
    if iWave < nWaves
        for a=1:3, for b=1:3
            h_base(:,:,:,a,b) = h_base(:,:,:,a,b) + H(a,b)*expPhase;
        end, end
    else
        for a=1:3, for b=1:3
            h_nth(:,:,:,a,b) = H(a,b)*expPhase;
        end, end
    end
end
fprintf("Done generating waves.\n\n");
%% Allocate frame storage
Cframes = cell(nFrames,1);
Lframes = cell(nFrames,1);
%% Compute each frame
fprintf("Computing frames...\n");
usePageMt = exist('pagemtimes','file')==2;
for f = 1:nFrames
    alpha = ampList(f);
    fprintf("  Frame %d / %d  (alpha = %.2f)\n", f, nFrames, alpha);
    % Field for current alpha
    h = h_base + alpha*h_nth;
    %% ---- Compute C-lines (Tr(h^2)=0) ----
    segmentsC = computeGWClines(h, X, Y, Z, x, y, z, usePageMt);
    Cframes{f} = segmentsC;
    %% ---- Compute L-points (n_h=0) ----
    Lpts = computeGWLpoints(h, X, Y, Z, x, y, z, usePageMt);
    Lframes{f} = Lpts;
end
fprintf("\nAll frames computed.\n");
%% Save animation data
save('GW_animation_data.mat', 'Cframes','Lframes','ampList','x','y','z','domainMin','domainMax');
fprintf("Saved to GW_animation_data.mat.\n");
function segmentsC = computeGWClines(h, X, Y, Z, x, y, z, usePageMt)
gridN = size(h,1);
%% Compute C = Tr(h^2)
if usePageMt
    hR = reshape(h,[],3,3);
    h2 = pagemtimes(hR,hR);
    Cvec = squeeze(h2(:,1,1) + h2(:,2,2) + h2(:,3,3));
    C = reshape(Cvec,gridN,gridN,gridN);
else
    C = zeros(gridN,gridN,gridN);
    for i=1:gridN, for j=1:gridN, for k=1:gridN
        Hloc = squeeze(h(i,j,k,:,:));
        C(i,j,k) = trace(Hloc*Hloc);
    end, end, end
end
ReC = real(C);
ImC = imag(C);
fvRe = isosurface(X,Y,Z,ReC,0);
segmentsC = [];
if isempty(fvRe.vertices)
    return;
end
FIm = griddedInterpolant({x,y,z},ImC,'linear','none');
V = fvRe.vertices;
vals = FIm(V(:,1),V(:,2),V(:,3));
faces = fvRe.faces;
for f=1:size(faces,1)
    idx = faces(f,:);
    p = V(idx,:);
    v = vals(idx);
    edgePts=zeros(3,3); nEdge=0;
    for e=1:3
        i1=e;
        i2=mod(e,3)+1;
        if v(i1)*v(i2)<0
            t = v(i1)/(v(i1)-v(i2));
            nEdge=nEdge+1;
            edgePts(nEdge,:) = p(i1,:) + t*(p(i2,:)-p(i1,:));
        end
    end
    if nEdge==2
        segmentsC = [segmentsC; edgePts(1,:); edgePts(2,:); NaN NaN NaN];
    end
end
end
function Lpoints = computeGWLpoints(h, X, Y, Z, x, y, z, usePageMt)
gridN = size(h,1);
Nx = zeros(gridN,gridN,gridN);
Ny = zeros(gridN,gridN,gridN);
Nz = zeros(gridN,gridN,gridN);
%% Compute n_h
if usePageMt
    hR = reshape(h,[],3,3);
    conjH = conj(hR);
    hT = permute(hR,[1 3 2]);
    S = pagemtimes(conjH,hT);
    n1 = 0.5*imag(S(:,2,3)-S(:,3,2));
    n2 = 0.5*imag(S(:,3,1)-S(:,1,3));
    n3 = 0.5*imag(S(:,1,2)-S(:,2,1));
    Nx = reshape(real(n1),gridN,gridN,gridN);
    Ny = reshape(real(n2),gridN,gridN,gridN);
    Nz = reshape(real(n3),gridN,gridN,gridN);
else
    for i=1:gridN, for j=1:gridN, for k=1:gridN
        Hloc = squeeze(h(i,j,k,:,:));
        S = conj(Hloc)*Hloc.';
        Nx(i,j,k) = 0.5*imag(S(2,3)-S(3,2));
        Ny(i,j,k) = 0.5*imag(S(3,1)-S(1,3));
        Nz(i,j,k) = 0.5*imag(S(1,2)-S(2,1));
    end, end, end
end
FNx = griddedInterpolant({x,y,z},Nx,'linear','nearest');
FNy = griddedInterpolant({x,y,z},Ny,'linear','nearest');
FNz = griddedInterpolant({x,y,z},Nz,'linear','nearest');
Nh2 = Nx.^2 + Ny.^2 + Nz.^2;
Nh_rms = sqrt(mean(Nh2(:)));
valTol = 1e-3 * Nh_rms;
Lpoints=[];
for ix=1:gridN-1
for iy=1:gridN-1
for iz=1:gridN-1
    idx = [
        ix   iy   iz;
        ix+1 iy   iz;
        ix   iy+1 iz;
        ix+1 iy+1 iz;
        ix   iy   iz+1;
        ix+1 iy   iz+1;
        ix   iy+1 iz+1;
        ix+1 iy+1 iz+1];
    cx = idx(:,1); cy=idx(:,2); cz=idx(:,3);
    NxC = Nx(sub2ind(size(Nx),cx,cy,cz));
    NyC = Ny(sub2ind(size(Ny),cx,cy,cz));
    NzC = Nz(sub2ind(size(Nz),cx,cy,cz));
    if ~(min(NxC)<0 && max(NxC)>0), continue; end
    if ~(min(NyC)<0 && max(NyC)>0), continue; end
    if ~(min(NzC)<0 && max(NzC)>0), continue; end
    x0 = mean(x([ix ix+1]));
    y0 = mean(y([iy iy+1]));
    z0 = mean(z([iz iz+1]));
    fNh2=@(r) (FNx(r(1),r(2),r(3))).^2 + ...
              (FNy(r(1),r(2),r(3))).^2 + ...
              (FNz(r(1),r(2),r(3))).^2;
    [rStar,fval] = fminsearch(fNh2,[x0,y0,z0],optimset('display','off'));
    if fval < valTol^2
        Lpoints=[Lpoints; rStar];
    end
end
end
end
%% Deduplicate nearby L-points
if ~isempty(Lpoints)
    d = (z(end)-z(1))/(gridN-1);
    keep = true(size(Lpoints,1),1);
    for i=1:size(Lpoints,1)
        if ~keep(i), continue; end
        d2=sum((Lpoints(i+1:end,:)-Lpoints(i,:)).^2,2);
        keep(i+find(d2<d^2))=false;
    end
    Lpoints=Lpoints(keep,:);
end
end

 
