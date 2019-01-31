function [OP,MM] = LvN_IPDG2D(potential_factor)

% Purpose: Set up the discrete Poisson matrix directly
%          using LDG. The operator is set up in the weak form
Globals2D;

% build local face matrices
massEdge = zeros(Np,Np,Nfaces);
Fm = Fmask(:,1); faceR = r(Fm); 
V1D = Vandermonde1D(N, faceR);  
massEdge(Fm,Fm,1) = inv(V1D*V1D');
Fm = Fmask(:,2); faceR = r(Fm); 
V1D = Vandermonde1D(N, faceR);  
massEdge(Fm,Fm,2) = inv(V1D*V1D');
Fm = Fmask(:,3); faceS = s(Fm); 
V1D = Vandermonde1D(N, faceS);  
massEdge(Fm,Fm,3) = inv(V1D*V1D');

% build local volume mass matrix
MassMatrix = invV'*invV;

% build DG derivative matrices
MM  = zeros(K*Np*Np, 3);  OP = zeros(K*Np*Np*(1+Nfaces), 3);  

% global node numbering
entries = (1:Np*Np)'; entriesMM = (1:Np*Np)'; 
for k1=1:K 
  if(~mod(k1,1000)) k1, end;
  rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';

  % Build local operators  
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;   Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;

  OP11 = J(1,k1)*(Dy'*MassMatrix*Dx/2 + Dx'*MassMatrix*Dy/2);

  % Build element-to-element parts of operator
  for f1=1:Nfaces
    k2 = EToE(k1,f1); f2 = EToF(k1,f1); 

    rows2 = ((k2-1)*Np+1:k2*Np)'*ones(1,Np); cols2 = rows2';
    
    fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + (1:Nfp);
    vidM = vmapM(fidM); Fm1 = mod(vidM-1,Np)+1;
    vidP = vmapP(fidM); Fm2 = mod(vidP-1,Np)+1;
    
    id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;
    % normals for id
    lnx = nx(id);  lny = ny(id); lsJ = sJ(id);
    hinv = max(Fscale(id), Fscale(1+(f2-1)*Nfp, k2));    

    Dx2 = rx(1,k2)*Dr + sx(1,k2)*Ds;   Dy2 = ry(1,k2)*Dr + sy(1,k2)*Ds;
    
    Dn1 = lnx*Dx  + lny*Dy ;
    Dn2 = lnx*Dx2 + lny*Dy2;

    mmE = lsJ*massEdge(:,:,f1);

    gtau = 90*2*(N+1)*(N+1)*hinv; % set penalty scaling
    switch(BCType(k1,f1))
      case {Dirichlet}
        OP11 = OP11 + ( gtau*mmE - mmE*Dn1 - Dn1'*mmE ); % ok
      case {Neuman}
        % nada 
      otherwise
        % interior face variational terms
        OP11        = OP11 + 0.5*( gtau*mmE - mmE*Dn1 - Dn1'*mmE );

        % off diagonal elements connecting the different elements k1 and k2
        % sharing the face f1  ---- > FLUX goes in here
        OP12 = zeros(Np);
        % set all rows in columns Fm2, with #Fm2 = Nfp (=6)
        OP12(:,Fm2) =             - 0.5*( gtau*mmE(:,Fm1) );    % nodes along f1 in element 1, "minus values"
        % set all columns in rows Fm1, with #Fm1 = Nfp (=6) and add
        % previously set values
        OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:) );
        OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*mmE(:, Fm1) );
        OP(entries(:), :) = [rows1(:), cols2(:), OP12(:)];  % so eine Art mesh-grid durch rows und cols erzeugt
        entries = entries + Np*Np;
    end 
  end      
  OP(entries(:), :)   = [rows1(:), cols1(:), OP11(:)];
  % apply the f(r,q) term, which is not the identity - care about
  % commutation with MassMatrix must be taken
  mM_times_factor = MassMatrix * potential_factor(rows1(:,1), cols1(1,:));
  MM(entriesMM(:), :) = [rows1(:), cols1(:), J(1,k1)*mM_times_factor(:)];
 % MM(entriesMM(:), :) = [rows1(:), cols1(:), J(1,k1)*MassMatrix(:)];
  entries = entries + Np*Np; entriesMM = entriesMM + Np*Np;
end  

OP   =   OP(1:max(entries)  -Np*Np,:);  OP   = myspconvert(OP, Np*K, Np*K, 1e-15);
MM   =   MM(1:max(entriesMM)-Np*Np,:);  MM   = myspconvert(MM, Np*K, Np*K, 1e-15);
return
