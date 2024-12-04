function [P1] = transition_probability(A, omega, T)
    a = -20;
    b = +20;
    L = b-a; 
    N = 10^3; 
    X1 = a+L*(0:N-1)/N;         %N grid points for FFT
    X2 = a+L*(1:N-1)/N;         %N-1 for DVR
    P = (2*pi/L)*[0:N/2-1,-N/2:-1];
    M = 10^3;
    dt = T/M;
    dx = L/N;
    
    UV_base = exp(-1i*(X1.^2/2)*dt/2);    % One-step propagator in position space
    UT = exp(-1i*(P.^2/2)*dt); 		% One-step propagator in momentum space
    
    % ground & first state
    
    % theory
    % n=0;
    % Hn = hermiteH(n, X1);
    % psi00 = (1 / sqrt(2^n * factorial(n) * sqrt(pi))) * exp(-X1.^2 / 2) .* Hn;
    %plot (X1(1:N),abs(psi00(1:N)).^2,'LineWidth',3,'LineStyle',':')
    %hold on
    
    %DVR method
    Tmatrix = zeros(N-1);       % Kinetic energy in DVR, (N-1)*(N-1) dimensional
    Vmatrix = (1/2)*diag(X2(1:N-1).^2);   % Potential energy in DVR
    nsum = 1:N-1;
    nsquare = (1:N-1).^2;
    c = (1/2)*(pi/L)^2*(2/N);
    for i1 = 1:N-1
        for i2 = i1:N-1
            Tmatrix(i1,i2) = c*sum(nsquare.*sin(nsum*pi*i1/N).*sin(nsum*pi*i2/N));
            Tmatrix(i2,i1) = Tmatrix(i1,i2);    % T is Hermitian
        end
    end
    Hamiltonian = Tmatrix+Vmatrix;            % Total Hamiltonian
    [vec,~] = eig(Hamiltonian);	% Eigenstates and eigenvalues
    vec0 = vec(1:N-1,1)./sqrt(sum(abs(vec(1:N-1,1)).^2*dx));  %ground state
    vec1 = vec(1:N-1,2)./sqrt(sum(abs(vec(1:N-1,2)).^2*dx));  %first excited state
    psi0 = [0;vec0]; %extend to N dimension
    psi1 = [0;vec1];
    % plot (X1(1:N),abs(psi0(1:N)).^2,'LineWidth',3,'LineStyle','--')
    % hold on
    
    %time evolution
    psi_tem=psi0'; %transpose if DVR
    %psi_tem=psi00
    P1_p=0;
    for m = 1:M
        t = m*dt;
        V_t = A * sin(X1).* cos(omega * t);
    
        %transition probability (1st perturbation)
        p_1=sum(psi1(1:N)'*diag(V_t(1:N))*psi0(1:N).*(L/N));
        pt=p_1*exp(1i*t)*dt;
        P1_p= P1_p+pt;
        %with RWA
        %V_rwa=A*sin(X1)
        %p_1rwa=sum(phi1(1:N)'*diag(V_rwa(1:N))*phi0(1:N).*L/N);
        %P1_rwa=p_1rwa*(exp(1i*(omega-1))-1)/(omega-1)
        UV_t = UV_base .* exp(-1i * V_t * dt / 2);
        psi_1 = UV_t.*psi_tem;
        phi_2 = fft(psi_1);
        phi_3 = UT.*phi_2;
        psi_3 = ifft(phi_3);
        psi_4 = UV_t.*psi_3;
        psi_tem = psi_4; 
    end
    
    psi=psi_tem(1:N)';       %final state
    P1 = [abs(psi1(1:N)'*psi(1:N)*L/N)^2,abs(P1_p)^2];    %transition probability
  
    %P1_perturbation_rwa=abs(P1_rwa)^2
    
    % plot (X1(1:N),abs(psi(1:N)).^2,'LineWidth',2)
    % hold on
end




