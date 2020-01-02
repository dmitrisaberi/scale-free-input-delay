% Function: protocol_design. 
%
function [epsilon, rho, taubar_max, P, K] = protocol_design(A,B,C,tau_bar)
    
    % get dimensions
    sizes = size(A);
    n = sizes(1);

    % get omega_max from the eigenvalues of A
    eigenvalues = eig(A);
    eigenvalues = eigenvalues';
    omega_max = -realmax;
    for lambda = eigenvalues
        if abs(real(lambda)) < 10^(-10)
            if imag(lambda) > omega_max
                omega_max = imag(lambda);
            end
        end
    end

    % if A is Hurwitz, set omega_max = 0
    if omega_max == -realmax
        omega_max = 0;
    end

    
    taubar_max = pi/(2*omega_max);
    
    if tau_bar >= taubar_max
        error("tau_bar is too large. Try a smaller value")
    end
    
    % fix rho
    rho = 2*(1/cos(tau_bar*omega_max));
    
    if rho <= 1/cos(tau_bar*omega_max)
        error("rho is too small. Try a larger value")
    end
    
    theta = (acos(1/rho)/tau_bar) - omega_max
    omega = omega_max + theta;
    omega_bar = max(norm(A) + 1, omega);
    length = omega_bar - omega;
    step_size = length/100;
    mu = sqrt(min(eig((A' + j*omega*eye(n))*(A - j*omega*eye(n)))));
    while abs(omega - omega_bar) > 10^(-5)*length
        for omega = omega_max + theta:step_size:omega_bar
            sv = sqrt(min(eig((A' + j*omega*eye(n))*(A - j*omega*eye(n)))));
            if mu >= sv
                mu = mu*0.999;
                break
            end
        end
    end
   
    % algorithm to compute epsilon
    % initialize values
    epsilon = 1;
    Q = eye(n);
    P = icare(A,[],Q,[],[],[],-B*B');
    value = norm(rho*B*B'*P);
    % repeat until value <= 0.99*mu
    while value > mu*0.99
        % update epsilon, value
        epsilon = 0.99*epsilon;
        Q = epsilon*eye(3);
        P = icare(A,[],Q,[],[],[],-B*B');
        value = norm(rho*B*B'*P);
    end
    
    % algorithm to get K. Change this if you want, but it will not
    % affect the solution too much.
    poles = [-1;-2;-3];
    K = place(A',C',poles);
    K = K';
    
end