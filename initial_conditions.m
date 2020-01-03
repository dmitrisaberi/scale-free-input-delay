    % Function: initial_conditions. This function serves to define the 
    % initial conditions for our time-delayed linear system. Since it is
    % infinite-dimensional (due to the delays), the initial conditions
    % are that of a function on [-tau_bar, 0]. We set it to zero
    % for all times prior to -tau_bar. You must give initial conditions
    % in the following format for our solver to work.
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
    %
    % The REAL initial conditions are ONLY for the agents x_1, ..., x_N.
    % The other ones (chi, xhat, x_r) do not affect the synchronization in 
    % any way. We suggest setting them to constant functions on
    % [-tau_bar, 0], for the sake of simplicity, as is done below.
    
     function s = initial_conditions(t)
        % define your parameters here *explicitly*. To get them,
        % run the above get_params function on your agent model
        % and delays.
        N = 3;
        n = 3;
        tau_bar = 0.3;
        % this parameter is 2 in the full-state case, 3 in the 
        % partial-state case. Delete it if you do not like this
        % syntax.
        alpha = 2;
         s = zeros(alpha*N*n+n,1);
           if t >= - tau_bar
               % Do not change the structure of the function.
               % *Only* edit inside of here, as it must be 0
               % outside of [-tau_bar, 0]. This DC signal initial 
               % condition is included here as a placeholder.
               for i = 1:alpha*N*n+n
                   s(i) = i;
               end
           end
     end
     
           