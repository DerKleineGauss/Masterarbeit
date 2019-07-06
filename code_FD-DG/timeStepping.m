function [v] = timeStepping(params, v)


rk4a = params.rk4a; rk4b = params.rk4b; rk4c = params.rk4c;
time = 0;

% Runge-Kutta residual storage
resv = zeros(params.Np*params.K* params.Ny, 1);

% get time-independet parts of systemmatrix
[L_glob, rhs, R, maxSpeed] = LVN_systemmatrix_timeIndependent(params);

% compute time step size
CFL=0.75; dt   = CFL/(maxSpeed*params.hy)*params.hx; dt = .5*dt;
Nsteps = ceil(params.FinalTime/dt); dt = params.FinalTime/Nsteps;

% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        
        % set potential
        B= functionB(params, timelocal);
        if(params.testing)
            figure(4)
            plot3(params.x_interface(:),params.y_interface(:), real(B(:)),'x')
            figure(5)
            plot3(params.x_interface(:),params.y_interface(:), imag(B(:)),'x')
        end
%         [G_glob] = LVN_systemmatrix_timeDependent(params, B, R);
        rhsv = (-L_glob-LVN_systemmatrix_timeDependent(params, B, R))*v + rhs;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        v = v+rk4b(INTRK)*resv;
    end
    % Increment time
    time = time+dt;
%     if (~mod(tstep,25))
%         plot_solution(params, v);
%     end
%     set ( line_handle, 'Xdata', reshape (x(2:Np ,:), [(Np-1)*K,1]), 'Ydata', reshape (u(2:Np ,:), [(Np-1)*K,1]));
end

end