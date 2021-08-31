

%%  Time vs state
% learned Q
load learned_Q Q;

thetas = -0.5*pi:1:pi;
qs = -1:0.1:1;
thetas_size = size(thetas,2);
qs_size = size(qs,2);

for theta = 1:thetas_size
    for q = 1:qs_size
    % initialize s
    x = [thetas(theta),qs(q)];%initial_state();

    % Initialize a new trial
    s = discretize_state(x, par);
    a = execute_policy(Q, s, par);

    % simulation time
    par.simtime = 4;     % Trial length PLACEHOLDER

    % bookkeeping (for ploting only)
    xa_history = zeros(ceil(par.simtime/par.simstep),numel(x));
    ta_history = 1:1:ceil(par.simtime/par.simstep);
    u_history  = zeros(ceil(par.simtime/par.simstep),numel(x));

    par.epsilon = 0; % EXPERIMENT
    % Inner loop: simulation steps
    for tt = 1:ceil(par.simtime/par.simstep)
        % Simulate a time step
        u = take_action(a, par);
        x = environment(x,u,par);

        % Save trace
        xa_history(tt, :) = x;
        u_history(tt, :) = u;
        s = discretize_state(x, par);
        a = execute_policy(Q, s, par);

        % Stop trial if state is terminal
        if is_terminal(s, par)
            disp('Goal reached');
            break
        end
    end
    end
end


function par = get_parameters(par)
    par.epsilon = 0.1;         % Random action rate
    par.gamma = 0.99;          % Discount rate
    par.alpha = 0.25;          % Learning rate
    par.theta_states = 100;    % angle discretization
    par.q_states = 100;        % angular velocity discretization
    par.actions = 15;          % Action discretization
    par.maxdeflection = 0.5;   % max fin deflection 
    par.trials = 5;            % Learning trials
end

function s = discretize_state(x, par)   
    % DISCRETIZED ANGULAR POSITION
    % Wrapped input between [0,2pi]
    wrap = abs(wrapTo2Pi(x(1)));
    % State in between [1,par.pos_states]
    if wrap >6.2831 % clip because at 2pi position was not discretized correctly
        theta = par.theta_states;
    else
        theta = floor((wrap *(par.theta_states/2)/pi)+1);
    end
    
    %DISCRETIZED ANGULAR VELOCITY
    %Clip input between [-pi,pi]
    if x(2)>=pi
        vwrap = pi;
    elseif x(2)<=-pi
        vwrap = -pi;
    else
        vwrap = x(2);
    end
    %State in between [1,par.vel_states]
    q = floor(interp1([-pi,pi],[1,par.q_states],vwrap));
    
    s = [theta,q];
end



