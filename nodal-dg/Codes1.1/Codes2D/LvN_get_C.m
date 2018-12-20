function [C] = get_C(r, N_q, q_interface, Delta_q)  
    C= zeros(N_q);
    for k=2:N_q-1
        C(k,k) = Delta_q/4 * (B(r, q_interface(k+1))+ B(r, q_interface(k))) ;
        C(k,k-1) = Delta_q/4 * (B(r, q_interface(k)));
        C(k,k+1) = Delta_q/4 * (B(r, q_interface(k+1))); 
    end
    C(1,1) = Delta_q/4 * (B(r, q_interface(1+1))+ B(r, q_interface(1))) ;
    C(1,2) = Delta_q/4 * (B(r, q_interface(1+1))); 
    C(N_q,N_q-1) = Delta_q/4 * (B(r, q_interface(N_q)));
    C(N_q,N_q) = Delta_q/4 * (B(r, q_interface(N_q+1))+ B(r, q_interface(N_q))) ;
end