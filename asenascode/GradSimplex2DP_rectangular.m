function [dmodedr, dmodeds] = GradSimplex2DP_rectangular(r,s,id,jd)

% function [dmodedr, dmodeds] = GradSimplex2DP_rectangular(r,s,id,jd)
% Purpose: Return the derivatives of the modal basis (id,jd) on the 2D simplex at (r,s).
fr = JacobiP(r, 0, 0, id);     dfr = GradJacobiP(r, 0, 0, id);
gs = JacobiP(s, 0, 0, jd);     dgs = GradJacobiP(s, 0, 0, jd);

dmodedr = dfr.*gs;
dmodeds = dgs.*fr;

return;
