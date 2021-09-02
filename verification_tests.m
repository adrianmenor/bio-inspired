%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots different parameters for code verification and data
% inspection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

test = 0; % 0 = system response, 1 = discretization

if test == 0
    %  open loop output
    x = [0.5,-0.9];
    par.simtime = 10;
    par.simstep = 0.05;
    x_history = zeros(par.simtime/par.simstep,2);
    t_history = zeros(par.simtime/par.simstep,1);
    u = 0;
    for t = 1:ceil(par.simtime/par.simstep)
        x = environment(x,u,par);
        t_history(t,:) = t;
        x_history(t,:) = x;
    end

    plot(t_history,x_history);
    
elseif test == 1
    %  Discretized state output
    par.theta_states = 100;
    par.q_states = 100;

    par.simstep = 0.05;
    
    % continuous state
    cs = [6.29,0];
    
    %discrete action
    fprintf('Discrete actions:');
    sP = discretize_state(cs,par)
    
    % obtained reward
    a = 1;
    fprintf('Reward: %d.\n',observe_reward(a, sP, par));
    
end

%% Functions being tested
    function s = discretize_state(x, par)   
        % DISCRETIZED ANGULAR POSITION
        % Wrap input between [0,2pi]
        wrap = abs(wrapTo2Pi(x(1)));
        % State in between [1,par.pos_states]
        if wrap >6.2831 % clip because at 2pi position was not discretized correctly
            theta = 31;
        else
            theta = floor((wrap *(par.theta_states/2)/pi)+1);
        end

        %DISCRETIZED ANGULAR VELOCITY
        %Clip input between [-5pi,5pi]
        if x(2)>=5*pi
            vwrap = 5*pi;
        elseif x(2)<=-5*pi
            vwrap = -5*pi;
        else
            vwrap = x(2);
        end
        %State in between [1,par.vel_states]
        q = floor(interp1([-5*pi,5*pi],[1,par.q_states],vwrap));

        s = [theta,q];
    end
    
function r = observe_reward(a, sP, par)
    % Calculates the reward for taking action a,
    % resulting in state sP.
    
    if sP(1)== 1%ceil(par.theta_states/2)% angle criteria
        if sP(2)==ceil(par.q_states/2)% q criteria
            r = 10;%reward
        else 
            r=0;%no reward
        end
    else
        r = 0;%no reward
    end
   
end














