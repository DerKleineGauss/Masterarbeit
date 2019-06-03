function [u] = Advec1D(u, FinalTime, params)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

time = 0;

% Runge-Kutta residual storage  
resu = zeros(params.Np,params.K); 

% advection speed
a = 2*pi;
a = a*10;

% compute time step size
xmin = min(abs(params.x(1,:)-params.x(2,:)));
CFL=0.75; dt   = CFL/(a)*xmin; dt = .5*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 

line_handle = plot(reshape(u,[params.K*params.Np, 1]));
axis([0 2 -1 1]);
% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + params.rk4c(INTRK)*dt;
        [rhsu] = AdvecRHS1D(u, timelocal, a, params);
        resu = params.rk4a(INTRK)*resu + dt*rhsu;
        u = u+params.rk4b(INTRK)*resu;
    end;
    % Increment time
    time = time+dt;
    set ( line_handle, 'Xdata', reshape (params.x(2:params.Np ,:), [(params.Np-1)*params.K,1]), 'Ydata', reshape (u(2:params.Np ,:), [(params.Np-1)*params.K,1]));
end;
return
