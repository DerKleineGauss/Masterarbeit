function [] = plotSolution_FVDG(u, params, realPart, nr_x, nr_y)

Lx = params.Lr_scaled;
Ly = params.Lq_scaled;
Kx = params.K;

y = params.y;
x = params.x;
y_interface = params.y_interface;

Np = params.Np;
Npy = params.Npy;   % number of nodes
Ny = params.Ny;       % equals 'N_q' in Lukas Paper : the number of cells

systemsize = K*Np*Ny;
hy = params.hy;

xlin= linspace(-Lx/2,Lx/2,nr_x);
ylin= linspace(-Ly/2,Ly/2,nr_y);
[X,Y]= meshgrid(xlin,ylin);
rhoint= zeros(nr_y,nr_x);

hx= abs(VX(2)-VX(1));
hy= abs(VY(2)-VY(1));



cellx= int64((xlin+Lx/2-mod(xlin+Lx/2,hx))/hx+1);
cellx(end)= cellx(end)-1;
celly= int64((ylin+Ly/2-mod(ylin+Ly/2,hy))/hy+1);
celly(end)= celly(end)-1;

offx= 0;
for j=1:Kx
    indx= find(cellx==j);
    offy=0;
    for k=1:Ny
        indy = find(celly==k);
        if ~isempty(indx) && ~isempty(indy)
            xloc= xlin(indx);
            yloc= ylin(indy);
            r= 2*mod(xloc+Lx/2,hx)/hx-1;
            s= 2*mod(yloc+Ly/2,hy)/hy-1;
            rhoc= rho((offy+1):(offy+Npy), (offx+1):(offx+Npx));
            Vsy= Vandermonde1D(Ny,s(:));
            Vsx= Vandermonde1D(Nx,r(:));
            rhoint(indy, indx)= Vsy*invVy*rhoc*invVx'*Vsx';
        end
        offy= offy+Ny;
    end
    offx= offx+Npx;
end

mesh(X,Y,real(rhoint))

end
