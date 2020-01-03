    % Function: initial_conditions. This function serves to define the 
    % initial conditions for our time-delayed linear system. Since it is
    % infinite-dimensional (due to the delays), the initial conditions
    % are that of a function on [-tau_bar, 0]. We set it to zero
    % for all times prior to -tau_bar. You must give initial conditions
    % in the following format for our solver to work. Make sure to define
    % them with the proper dimensions, depending on N,n and if it is
    % full or partial-state coupling.
    %
    % Full State: You will define a function from 
    % [- tau_bar, 0] -> R^{2*N*n + n}. The first N*n variables of the
    % output are the initial conditions of the agents, x_1, ..., x_N.
    % The next N*n variables are that of chi_1, ..., chi_N, 
    % whose initial conditions can be set to whatever you like, 
    % so long as they are continuous. There are no restrictions
    % on the choice of initial conditions for the final n variables, x_r.
    %
    % Partial State: You will define a function from 
    % [-tau_bar, 0] -> R^{3*N*n + n}. Similar to above, the first
    % N*n variables are x_1, ..., x_N, which is where your initial 
    % conditions go. The next 2*N*n variables, xhat_1, ..., xhat_N, 
    % chi_1, ..., chi_N can have any arbitrary continuous initial
    % conditions. The final n variables, x_r, again have no restrictions.
  
    function s = initial_conditions(t)
        % define your parameters here *explicitly*. tau_bar
        % is always necessary to define (since the function
        % is only nonzero on [-tau_bar,0]), but N,n are only
        % necessary if you refer to them below.
        N = 3;
        n = 3;
        tau_bar = 0.3;
        % this parameter is 2 in the full-state case, 3 in the 
        % partial-state case. Delete it if you do not make
        % use of it below.
        alpha = 2;
        s = zeros(alpha*N*n+n,1);
        if t >= - tau_bar
               % Do not change the structure of the function.
               % Insert your initial conditions *here*, in this if
               % statement, and don't change the structure
               % of the function (including the initialization
               % of s above).
            for i = 1:alpha*N*n+n
                   s(i) = i;
            end
        end
    end
     
           
