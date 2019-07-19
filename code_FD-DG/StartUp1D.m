function [out_params] = StartUp1D(params)
    % Purpose : Setup script, building operators, grid, metric and connectivity for 1D solver.     

    % Definition of constants

    params.NODETOL = 1e-10;
    params.Np = params.N+1; 
    params.Nfp = 1; 
    params.Nfaces = 2;

    % Compute basic Legendre Gauss Lobatto grid
    params.r = JacobiGL(0,0,params.N);

    % Build reference element matrices
    params.V  = Vandermonde1D(params, params.r); params.invV = inv(params.V);
    params.MassMatrix = params.invV'*params.invV;
    params.Dr = Dmatrix1D(params);

    % Create surface integral terms
    params.LIFT = Lift1D(params);

    % build coordinates of all the nodes
    va = params.EToV(:,1)'; vb = params.EToV(:,2)';
    params.x = ones(params.N+1,1)*params.VX(va) + 0.5*(params.r+1)*(params.VX(vb)-params.VX(va));

    % calculate geometric factors
    [params.rx,params.J] = GeometricFactors1D(params);
    if (params.N == 0)
        params.J = 0*params.x + params.hx/2;
    end

    % Compute masks for edge nodes
    fmask1 = find( abs(params.r+1) < params.NODETOL)'; 
    fmask2 = find( abs(params.r-1) < params.NODETOL)';
    params.Fmask  = [fmask1;fmask2]';
    if (params.N == 0)
        params.Fmask = 1;
    end
    params.Fx = params.x(params.Fmask(:), :);

    % Build surface normals and inverse metric at surface
    [params.nx] = Normals1D(params);
    params.Fscale = 1./(params.J(params.Fmask,:));

    % Build connectivity matrix
    [params.EToE, params.EToF] = Connect1D(params);

    % Build connectivity maps
    if (params.N > 0)
        [params] = BuildMaps1D(params);
    end
    
    for k1=1:params.K
      params.loc2glb(:,k1) = ((k1-1)*params.Np+1:k1*params.Np)';
    end
    
    out_params = params;
end
