function [OP] = PoissonIPDG2D_rectangular(params, B)

% Purpose: Set up the discrete Poisson matrix directly
%          using LDG. The operator is set up in the weak form

Nfaces = params.Nfaces;
Np = params.Np;
Nfp = params.Nfp;
Nx = params.Nx;
Ny = params.Ny;
K = params.K;
Fmask = params.Fmask;
vmapM = params.vmapM;
vmapP = params.vmapP;
Fscale = params.Fscale;
r = params.r;
s = params.s;
MassMatrix = params.MassMatrix;
Dr = params.Dr;
Ds = params.Ds;
rx = params.rx;
ry = params.ry;
sx = params.sx;
sy = params.sy;
nx = params.nx;
ny = params.ny;
J = params.J;
sJ = params.sJ;
EToE = params.EToE;
EToF = params.EToF;
BCType = params.BCType;
Dirichlet = params.Dirichlet;
Neuman = params.Neuman;

% build local face matrices
massEdge = zeros(Np,Np,Nfaces);
Fm = Fmask(:,1); faceR = r(Fm);
V1D = Vandermonde1D(Nx, faceR); massEdge(Fm,Fm,1) = inv(V1D*V1D');
Fm = Fmask(:,2); faces = s(Fm);
V1D = Vandermonde1D(Ny, faces); massEdge(Fm,Fm,2) = inv(V1D*V1D');
Fm = Fmask(:,3); facer = r(Fm);
V1D = Vandermonde1D(Nx, facer); massEdge(Fm,Fm,3) = inv(V1D*V1D');
Fm = Fmask(:,4); faceS = s(Fm);
V1D = Vandermonde1D(Ny, faceS); massEdge(Fm,Fm,4) = inv(V1D*V1D');


% build DG derivative matrices
% MM  = zeros(K*Np*Np, 3);
% OP = zeros(K*Np*Np*(1+Nfaces), 3);

% global node numbering
entries = (1:Np*Np)';
indices_row = zeros(K*Np*Np*(1+Nfaces),1);
indices_col = zeros(K*Np*Np*(1+Nfaces),1);
values = zeros(K*Np*Np*(1+Nfaces),1);

for k1=1:K
  if(~mod(k1,1000)) k1, end;
  rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np);
  cols1 = rows1';

  % Build local operators
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;
  Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;

  OP11 = J(1,k1)*(Dx'*MassMatrix*Dy + Dy'*MassMatrix*Dx+MassMatrix*B(rows1(:,1), cols1(1,:)));

  % Build element-to-element parts of operator
  for f1=1:Nfaces
    k2 = EToE(k1,f1);
    f2 = EToF(k1,f1);

    rows2 = ((k2-1)*Np+1:k2*Np)'*ones(1,Np); cols2 = rows2';

    fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + (1:Nfp);
    vidM = vmapM(fidM); Fm1 = mod(vidM-1,Np)+1;
    vidP = vmapP(fidM); Fm2 = mod(vidP-1,Np)+1;

    id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;
    lnx = nx(id);  lny = ny(id); lsJ = sJ(id);
    hinv = max(Fscale(id), Fscale(1+(f2-1)*Nfp, k2));

    Dx2 = rx(1,k2)*Dr + sx(1,k2)*Ds;   Dy2 = ry(1,k2)*Dr + sy(1,k2)*Ds;

    Dn1 = lnx*Dx  + lny*Dy;
    Dn2 = lnx*Dx2 + lny*Dy2;

    Dn3 = lny*Dx + lnx*Dy;

    mmE = lsJ*massEdge(:,:,f1);

    gtau = 100*2*(Nx+1)*(Nx+1)*hinv; % set penalty scaling
    %fprintf('BCType(k1,f1): %d\n',BCType(k1,f1))
    switch(BCType(k1,f1))
      case {Dirichlet}
	OP11 = OP11 + ( gtau*mmE - mmE*Dn1 - Dn3'*mmE ); % ok
      case {Neuman}
        % nada
      otherwise
	% interior face variational terms
	OP11        = OP11 + 0.5*( gtau*mmE - mmE*Dn1 - Dn3'*mmE );

	OP12 = zeros(Np);
	OP12(:,Fm2) =             - 0.5*( gtau*mmE(:,Fm1) );
	OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:) );
	OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn3'*mmE(:, Fm1) );
    indices_row(entries(:)) = rows1;
    indices_col(entries(:)) = cols2;
    values(entries(:)) = OP12;
	entries = entries + Np*Np;
    end
  end
  indices_row(entries(:)) = rows1;
  indices_col(entries(:)) = cols1;
  values(entries(:)) = OP11;
%   MM(entriesMM(:), :) = [rows1(:), cols1(:), J(1,k1)*MassMatrix(:)];
  entries = entries + Np*Np;
  %fprintf('k1: %d\tmin: %d\tmax: %d\n\n',k1,min(entries),max(entries))
end

indices_row=indices_row(1:max(entries)  -Np*Np,:);
indices_col=indices_col(1:max(entries)  -Np*Np,:);
values=values(1:max(entries)  -Np*Np,:);
OP = [indices_row(:), indices_col(:), values(:)];
OP   = myspconvert(OP, Np*K, Np*K, 1e-15);


% MM   =   MM(1:max(entriesMM)-Np*Np,:);  MM   = myspconvert(MM, Np*K, Np*K, 1e-15);

return
