
% Input Delay Solver. This function accepts the following parameters:

% 1. The agent model (A, B, C). Types: matrices.
% 2. The input delays (taus). Type: array.
% 3. The adjacency matrix of the graph (A_script). Type: matrix.
% 4. The set of leader agents (leader_set). Type: array.
% 5. The initial conditions (initial_conditions). Type: function.
% (Note: define these in a separate m-file, a template for which is given.
% the necessary form of the initial condition function is described 
% in the paper).
% 6. Tmax: The period of integration, i.e., the solution is given
% on the interval [0, Tmax].

% From this information, input_delay_solver implements the protocol 
% described in Saberi, et al. Specifically, the linear system is solved, 
% and synchronization of the agents is achieved. If the given C is the identity,
% full-state coupling is implemented. Otherwise, the protocol for 
% partial-state coupling is enacted. The core solver used
% is dde23, which chooses a mesh (i.e., a time series) and gets the 
% approximate state value at each time. In the end, input_delay_solver
% returns 4 objects:

% 1. t (the time series). Type: array.
% 2. x (the state). Type: matrix (specifically, an n by T matrix, where n
% is the dimension of the state and T is the number of elements in the time
% series. The columns of this matrix give the state values at the corresponding time step).
% 3. x_r (the state of the exosystem). Type: matrix (same structure as
% above).
% 4. u (the input). Type: matrix (same structure as above).
%
% Some things you can do after solving: 
% 1. plot the states: call plot(t,x).
% 2. plot xtilde (which should -> 0). Note x = [x_1; x_2; ...; x_N], so
% xtilde = [x_1 - x_r; x_2 - x_r; ...; x_N - x_r]. Then, call
% plot(t,xtilde).

function [t x x_r u] = input_delay_solver(A,B,C, taus, initial_conditions, ...
    A_script, leader_set, T_max)

    tau_bar = max(taus);
    tspan = [0 T_max];
    
    % extract dimensions from data (i.e., N, n, q)
    sz1 = size(A);
    N = sz1(1);
    sz2 = size(A_script);
    n = sz2(1);
    sz3 = size(C);
    q = sz3(1);

    % get leader from leader_set
    leader = zeros(N,1);
    for i = 1:N
        if ismember(i, leader_set)
            leader(i) = 1;
        end
    end
    
    % get Laplacian
    L = zeros(N,N);
    for i = 1:N 
        v = 0;
        for j = 1:N
            v = v + A_script(i,j);
        end
        L(i,i) = v;
    end
    
    for i = 1:N
        for j = 1:N
            if i ~= j
                L(i,j) = -A_script(i,j);
            end
        end
    end
    L_bar = L + diag(leader);
    
    % design your protocol, using the provided function
    [epsilon, rho, taubar_max, P, K] = protocol_design(A,B,C,tau_bar);

    % Dtilde = -rho * [zeros(N*n, N*n) kron(eye(N), B*B'*P); ...
        % zeros(N*n, N*n) kron(eye(N), B*B'*P)];
    % Atilde = [kron(eye(N), A) zeros(N*n, N*n); kron(L_bar, eye(n))...
    % kron(eye(N), A) - kron(L_bar, eye(n))];
    
    Dtilde = 0
    Atilde = 0

    % frame of algorithm
    if C == eye(n)
    % full-state coupling
        Ltilde = kron(leader,eye(n));
        Dtilde = -rho * [zeros(N*n, N*n) kron(eye(N), B*B'*P); ...
        zeros(N*n, N*n) kron(eye(N), B*B'*P)];
        Dtilde = [Dtilde zeros(2*N*n,n); zeros(n,2*N*n) zeros(n,n)];
        Atilde = [kron(eye(N), A) zeros(N*n, N*n) zeros(N*n, n); ...
            kron(L_bar, eye(n)) kron(eye(N), A) - kron(L_bar, eye(n)) ...
            -Ltilde; zeros(n, 2*N*n) A];
    else
    % partial-state coupling
        Atilde = [kron(eye(N), A) zeros(N*n,N*n) zeros(N*n,N*n) ...
            zeros(N*n,n); kron(L_bar,K*C) kron(eye(N), A - K*C) ...
            zeros(N*n,N*n) -kron(leader,K*C); zeros(N*n,N*n) ...
            eye(n*N) kron(eye(N),A) - kron(L_bar,eye(n)) zeros(n*N,n); ...
            zeros(n,N*n) zeros(n,N*n) zeros(n,N*n) A];
        Dtilde = -rho*[zeros(N*n,N*n) zeros(N*n,N*n) kron(eye(N), B*B'*P) ...
            zeros(N*n,n); zeros(N*n,N*n) zeros(N*n,N*n) kron(L_bar, B*B'*P) ...
            zeros(N*n,n); zeros(N*n,N*n) zeros(N*n,N*n) kron(eye(N), B*B'*P) ...
            zeros(N*n,n); zeros(n,N*n) zeros(n,N*n) zeros(n,N*n) ...
            zeros(n,n)];
    end

    % get Atilde
   % Leader_matrix = [zeros(N*n, N*n) zeros(N*n, N*n); zeros(N*n) ...
    % kron(diag(leader), eye(n))];
    % Atilde = Atilde - Leader_matrix;
    
    % Function: my_ddefun. This function serves to define
    % the time-delayed linear system we wish to simulate,
    % with full-state coupling. It accepts a time parameter t,
    % state variable x (with dimensions 2*N*n, as it represents the
    % column vector (x_1, ..., x_N, chi_1, ..., chi_N), with the x_i
    % and chi_i defined in the paper), and returns a 2*N*n-dimensional
    % vector dydt, which is the right-hand side of the linear system.
    function dydt = my_ddefun(t,x,Z)
        % full state
        if C == eye(n)
            dydt = zeros(2*N*n+n,1);
            rhs = zeros(2*N*n+n, 1);
        else % partial state
            dydt = zeros(3*N*n+n,1);
            rhs = zeros(3*N*n+n, 1);
        end
        rhs = Atilde*x;
        for j = 1:N
            R_j = blkdiag(zeros((j-1)*n, (j-1)*n), eye(n), ...
                zeros((N-1)*n, (N-1)*n), eye(n), zeros((N-j)*n, (N-j)*n));
            if C == eye(n) % full state
                R_j = [R_j zeros(2*N*n,n); zeros(n,2*N*n) zeros(n,n)];
            else % partial state
                R_j = [R_j zeros(2*N*n,n); zeros(n,2*N*n) zeros(n,n)];
                R_j = [zeros(N*n, 3*N*n +n); zeros(2*N*n + n, N*n) R_j];
            end
            
            A_j = Dtilde*R_j;
            rhs = rhs + A_j*Z(:,j);
        end
        dydt = rhs;
    end

    % solve the system
    sol = dde23(@my_ddefun, taus, @initial_conditions, tspan);
  
    t = sol.x;
    x = sol.y(1:N*n,:);
    chi = 0;
    if C == eye(n)
        % full state
        x_r = sol.y(2*N*n+1:2*N*n+n,:);
        chi = sol.y(N*n+1:2*N*n,:);
    else 
       %  partial state
        x_r = sol.y(3*N*n+1:3*N*n + n, :);
        chi = sol.y(2*N*n+1:3*N*n,:);
    end
    u = -rho*kron(eye(N), B*B'*P)*chi;

end
