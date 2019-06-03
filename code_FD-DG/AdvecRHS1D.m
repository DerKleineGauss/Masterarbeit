function [rhsu] = AdvecRHS1D(u,time, a, params)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

% form field differences at faces
alpha=1;
du = zeros(params.Nfp*params.Nfaces,params.K); 
du(:) = (u(params.vmapM)-u(params.vmapP)).*(a*params.nx(:)-(1-alpha)*abs(a*params.nx(:)))/2;

% impose boundary condition at x=0
uin = -sin(a*time);
du (params.mapI) = (u(params.vmapI)- uin ).*(a*params.nx(params.mapI)-(1-alpha)*abs(a*params.nx(params.mapI)))/2;
du (params.mapO) = 0;

% compute right hand sides of the semi-discrete PDE
rhsu = -a*params.rx.*(params.Dr*u) + params.LIFT*(params.Fscale.*(du));
return
