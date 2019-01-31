function varargout = getHartreeFockPotential(eps_V, V, n, n_cell_r, dr, Nd_V)

% Implementation of the numerical solution of the poisson equation

% get the constants

loadConstants;

% acceleration of convergence

V_ref = kB*T/q;

P = diag(-(eps_V(1:n_cell_r)+eps_V(2:n_cell_r+1))./(dr.^2)-q*n./V_ref, 0) + ...
    diag( eps_V(2:n_cell_r)./(dr.^2), -1) + ...
    diag( eps_V(2:n_cell_r)./(dr.^2), +1);

R =q*(n.*(1-V'./V_ref)-Nd_V).';

% R(1) = R(1) - V(1)*P(2,1);
% R(n_cell_r) = R(n_cell_r) - V(n_cell_r)*P(n_cell_r-1, n_cell_r);

R(1) = R(1);
R(n_cell_r) = R(n_cell_r) - V_apl*P(n_cell_r-1, n_cell_r);

varargout{1} = P\(R);