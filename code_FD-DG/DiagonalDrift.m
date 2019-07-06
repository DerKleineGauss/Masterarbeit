function [eigs_A, R] = DiagonalDrift(params)   
    % Decomposition A = R \Lambda R' -> u=Rv for A with zero dirichlet bc's
    Nq = params.Ny;
    Npy = Nq+1;
    a = 0;
    b = +1i/2;c = -1i/2; factor = sqrt(b*c);

    % sqrt(bc) = sqrt(b/c) = 1i
    eigs_A = zeros(Npy-1,1);
    R = zeros(Npy-1);
    norm_R=0;
    if(mod(Npy,2) == 0)
        error("Number of nodes (cells) must in y direction be odd (even).");
    end

    k=1:Nq;
    for j=1:Npy-1
        eigs_A(j) = (a - sign(imag(c)) * 2*factor*cos(pi*j/((Npy-1)+1)));
        R(:,j) = (b/c).^(k/2).*sin(k*pi*j/Npy);
        norm_R = norm_R+R(1,j)^2;
    end
    R = 1/sqrt(norm_R) * R;
    Lambda = diag(eigs_A);

    % Tests
    if (params.testing)
            A = full(gallery('tridiag',(Npy-1),b,0,c)) ;
        display(['|| id - R_invers R || = ', num2str(norm(eye(Npy-1)-R'*R))]);
        Test = (A - R*Lambda*R');
        display(['|| A - R \Lambda R_invers || = ' , num2str(norm(Test))]);
        display(['|| (A - \lambda_1 * id)*r_1 || = ', num2str(norm((A-eigs_A(1)*eye(Npy-1))*R(:,1)) )])
    end
        

    %% periodic boundary conditions:
%     a = +1i/2;
%     b = -1i/2;
%     Nq = params.Ny;
%     Npy = Nq+1;
%     normR = sqrt(Nq);
%     if(mod(Npy,2) == 0)
%         error("Number of nodes (cells) in y direction must be odd (even).");
%     end
%     eigs_A = zeros(Nq,1);
%     R = zeros(Nq);
%     
%     n = fftshift(0:Nq-1);
%     k=1:Nq;
%     for j=1:Nq
%         eigs_A(j) = -2i*a*sin(2*pi*n(j)/Nq);
%         R(:,j) = exp(2i*pi*n(j).*k./Nq)/normR;
%     end
% 
%         % Tests
%         if (params.testing)
%             A = full(gallery('tridiag',(Npy-1),a,0,b)) ;
%             A(1,Nq) = a;
%             A(Nq,1) = b;
%             Lambda = diag(eigs_A);
%             display(['|| id - R_invers R || = ', num2str(norm(eye(Npy-1)-R'*R))]);
%             Test = (A - R*Lambda*R');
%             display(['|| A - R \Lambda R^T || = ' , num2str(norm(Test))]);
%             display(['|| (A - \lambda_1 * id)*r_1 || = ', num2str(norm((A-eigs_A(1)*eye(Npy-1))*R(:,1)) )])
%         end