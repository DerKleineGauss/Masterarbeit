function [] = plotDensity_FVDG(v, params, R, n_linspace, fid, timelocal)

% function [] = plotDensity_FVDG(v, params, R, n_interpol)
% Purpose : Plot density function of finite volume / dg solution v which is
%           in eigenspace
%           R :         matrix of eigenvectors
%           n_linspace: size of linspace where v shall be evaluated
%           fid:        id of figure to get with fid =figure;

Lx = params.Lr_scaled;
Ly = params.Lq_scaled;
Kx = params.K;

y = params.y;
x = params.x;
y_interface = params.y_interface;

Np = params.Np;
Npy = params.Npy;   % number of nodes
Ny = params.Ny;       % equals 'N_q' in Lukas Paper : the number of cells

hy = params.hy;

v = reshape(v,params.Np*params.K, params.Ny);
u = strangeMatrixMultiplication(R , v);
u = reshape(u,params.Np*params.K, params.Ny);

% y = 0 is always the interface between cell Ny/2 and cell Ny/2+1
u_0 = (u(:,Ny/2)+u(:,Ny/2 + 1)) / 2;

xlin= linspace(-Lx/2,Lx/2,n_linspace);

rhoint= 0*xlin;

hx= params.VX(2)-params.VX(1);



cellx= int64((xlin+Lx/2-mod(xlin+Lx/2,hx))/hx+1);
% cellx(end)= cellx(end)-1;

offx= 0;
for j=1:Kx
    indx= find(cellx==j);
    if ~isempty(indx)
        xloc= xlin(indx);
        r= 2*mod(xloc+Lx/2,hx)/hx-1;
        uc= u_0((offx+1):(offx+Np));
        Vsx= Vandermonde1D(params, r(:));
        invV = params.invV;
        rhoint(indx)= uc'*invV'*Vsx';
    end

    offx= offx+Np;
end

figure(fid)
clf(fid)
yyaxis left
plot(xlin,real(rhoint),'.')
hold on
plot(xlin,real(rhoint))
yyaxis right
plot(xlin, Potential(xlin, timelocal, params))
hold off

end
