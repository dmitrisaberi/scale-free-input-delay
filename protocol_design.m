% Protocol Design. This function accepts the following parameters:
% 
% 1. The agent model: A, B, C. Types: matrices.
% 2. An upper bound on the delays: tau_bar. Type: positive number.
% 
% From this information, protocol_design picks parameters that will
% go into the protocol described in the Scale-Free input delay paper. 
% Our metric of performance is given by speed of synchronization, i.e., 
% how quickly the agents x_i converge to the exosystem x_r. Typically,
% larger epsilon lead to better performance. As such, we have designed
% our algorithm with an eye towards larger epsilon. The algorithm is
% designed as follows. First, we find the maximum upper bound on delays
% (taubar_max) through the solvability condition given in the paper
% (Equation (16)). Then, we check that the given tau_bar is strictly less
% than taubar_max. If the problem is solvable, we design the protocol,
% returning the following parameters:

% 1. epsilon. Type: positive number.
% 2. rho. Type: positive number.
% 3. taubar_max. Type: positive number.
% 4. P. Type: positive definite matrix.
% 5. K. Type: matrix.

% NOTE: This algorithm is meant to generate valid parameters (epsilon, rho) 
% given any agent model and valid upper bound on delays. We make no claim
% that this choice is optimal for convergence. In fact, you may often
% find that choosing epsilon larger than we give here will make performance
% much better. What we DO guarantee is that you need not choose epsilon any
% smaller than what this function returns. We welcome any improvements 
% upon this algorithm to get a larger allowable epsilon.
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

    % get taubar_max
    taubar_max = pi/(2*omega_max);
    
    if tau_bar >= taubar_max
        error("tau_bar is too large. Try a smaller value")
    end
    
    % fix rho
    rho = 2*(1/cos(tau_bar*omega_max));
    
    if rho <= 1/cos(tau_bar*omega_max)
        error("rho is too small. Try a larger value")
    end
    
    % Compute theta
    theta = (acos(1/rho)/tau_bar) - omega_max
    
    % algorithm to compute a non-conservative mu
    omega = omega_max + theta;
    omega_bar = max(norm(A) + 1, omega);
    length = omega_bar - omega;
    step_size = length/100;
    mu = sqrt(min(eig((A' + j*omega*eye(n))*(A - j*omega*eye(n)))));
    while abs(omega - omega_bar) > 10^(-5)*length
        % search the interval (omega_max + theta, omega_bar) to see
        % if the given mu minimizes sv
        for omega = omega_max + theta:step_size:omega_bar
            sv = sqrt(min(eig((A' + j*omega*eye(n))*(A - j*omega*eye(n)))));
            if mu >= sv
                % mu does not minimize sv; update mu
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