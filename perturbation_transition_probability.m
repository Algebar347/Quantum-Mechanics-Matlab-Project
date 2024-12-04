function [P1_perturbation] = perturbation_transition_probability(T, A, omega)
    % Parameters
    a = -20;
    b = +20;
    L = b-a;
    N = 10^3;
    X1 = a+L*(0:N-1)/N; %N grid points for FFT
    X2 = a+L*(1:N-1)/N; %N-1 for DVR
    P = (2*pi/L)*[0:N/2-1,-N/2:-1];
    
    % Time parameters
    M = 10^3;
    dt = T/M;
    dx = L/N;
    
    % Propagators
    UV_base = exp(-1i*(X1.^2/2)*dt/2); % One-step propagator in position space
    UT = exp(-1i*(P.^2/2)*dt); % One-step propagator in momentum space
    
    % DVR method for Hamiltonian
    Tmatrix = zeros(N-1);
    Vmatrix = (1/2)*diag(X2(1:N-1).^2);
    nsum = 1:N-1;
    nsquare = (1:N-1).^2;
    c = (1/2)*(pi/L)^2*(2/N);
    
    for i1 = 1:N-1
        for i2 = i1:N-1
            Tmatrix(i1,i2) = c*sum(nsquare.*sin(nsum*pi*i1/N).*sin(nsum*pi*i2/N));
            Tmatrix(i2,i1) = Tmatrix(i1,i2);
        end
    end
    
    Hamiltonian = Tmatrix+Vmatrix;
    [vec,~] = eig(Hamiltonian);
    
    % Ground and first excited states
    vec0 = vec(1:N-1,1)./sqrt(sum(abs(vec(1:N-1,1)).^2*dx));
    vec1 = vec(1:N-1,2)./sqrt(sum(abs(vec(1:N-1,2)).^2*dx));
    
    psi0 = [0;vec0];
    psi1 = [0;vec1];
    
    % Time evolution
    psi_tem = psi0';
    P1_p = 0;
    
    for m = 1:M
        t = m*dt;
        V_t = A * sin(X1).* cos(omega * t);
        
        % Perturbation calculation
        p_1 = sum(psi1(1:N).*V_t(1:N).*psi0(1:N).*L/N);
        pt = p_1*exp(1i*t)*dt;
        P1_p = P1_p+pt;
        
        % Time propagation
        UV_t = UV_base .* exp(-1i * V_t * dt / 2);
        psi_1 = UV_t.*psi_tem;
        phi_2 = fft(psi_1);
        phi_3 = UT.*phi_2;
        psi_3 = ifft(phi_3);
        psi_4 = UV_t.*psi_3;
        psi_tem = psi_4;
    end
    
    psi = psi_tem(1:N)';
    
    % Calculate transition probabilities
    P1_perturbation = abs(P1_p).^2;
end