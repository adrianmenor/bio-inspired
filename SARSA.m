%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code runs the Sarsa algorithm on a discrete state and action space,
% and outputs the learned Q table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Clearing workspace
clear all
close all
clc
rand('seed', 9); %7

%%  Sarsa
tic
par.simtime = 10;     % Trial length
par.simstep = 0.05;   % Simulation time step

% Obtain Sarsa parameters
par = get_parameters(par);

% Initialize Q(s,a) arbitrarily
% Q = init_Q(par);
load learned_Q Q; % PLACEHOLDER
Q = Q;            % PLACEHOLDER
% Initialize bookkeeping (for plotting only)
ra = zeros(par.trials, 1);
ra_cumulative = zeros(par.trials, 1);
tta = zeros(par.trials, 1);

% Outer loop: trials
for ii = 1:par.trials

    % Initializing the inner loop
    x = initial_state();%initialize s
    s = discretize_state(x, par);%discretize initial s
    a = execute_policy(Q, s, par);%choose a from s
    % Inner loop: simulation steps
    for tt = 1:ceil(par.simtime/par.simstep)

        % take action a
        u  = take_action(a, par);
        % Apply deflection and obtain new state
        % x  : state 
        % u  : action     
        x = environment(x,u,par);
        
        % learn:
        % using s for discretized state        
        sP = discretize_state(x, par);
        
        reward = observe_reward(u, sP, par);       
        aP = execute_policy(Q, sP, par); %choose a from sP

        Q = update_Q(Q, s, a, reward, sP, aP, par);
        s = sP; %update state
        a = aP; %updated action
        % Keep track of cumulative reward
        ra_cumulative(ii) = sum(ra)+reward;
        ra(ii) = ra(ii)+reward;
        
        % check termination condition
        if is_terminal(s, par)  
            fprintf('Goal reached in trial: %d.\n',ii);
            break
        end
    end

    tta(ii) = tta(ii) + tt*par.simstep;

end
toc

% saving learned Q value function
par.Q = Q;
save('learned_Q.mat','Q');

%%  Testing controller after learning
% learned Q
Q = par.Q;

% initialize s
x = [-pi,1];

% Initialize a new trial
s = discretize_state(x, par);
a = execute_policy(Q, s, par);

% simulation time
par.simtime = 12;     % Trial length PLACEHOLDER

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

%% Plotting 
% subplot(2,1,1);
% % plotting cumulative reward
% xaxis = 1:1:size(ra,1);
% plot(xaxis,ra_cumulative);
% title('Cumulative reward');
% xlabel('Trial [-]');
% ylabel('Reward [-]');
% 
% subplot(2,1,2); 
% plotting states over time
plot(ta_history,xa_history(:,1),'r','LineWidth',1.5);
hold on
plot(ta_history,xa_history(:,2),'b','LineWidth',1.5);
plot(ta_history,u_history,'k','LineWidth',.5);
hold off
title('State evolution and control input','FontSize',40);
ax = gca;
ax.FontSize = 30;
legend('\Theta [rad]','q [rad/s]','\delta_e [rad]','FontSize',30)
xlabel('Simulation steps [-]','FontSize',30);
ylabel('State [-]','FontSize',30);

%%  Functions
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

function x0 = initial_state()
    % random initial state
%     x0 = [2*rand(),-1+2*rand()]; 
    x0 = [-50*pi+100*pi*rand(),-50*pi+100*pi*rand()];
end

function Q = init_Q(par)
    %Initial Q = +1
    Q = ones(par.theta_states, par.q_states, par.actions);
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

function u = take_action(a, par)
    % Calculates the proper torque for action a. This cannot
    % exceed par.maxdeflection.
    u = (interp1([1,par.actions],[-par.maxdeflection,par.maxdeflection],a));
end

function r = observe_reward(a, sP, par)
    % Calculates the reward for taking action a,
    % resulting in state sP.
    
    if sP(1)== 1% angle criteria
        if sP(2)==ceil(par.q_states/2)% q criteria
            r = 10;%reward
        else 
            r=0;%no reward
        end
    else
        r = 0;%no reward
    end
   
end

function t = is_terminal(sP, par)
    %  Returns 1 if state sP is terminal, 0 otherwise.
    if sP(1)== 1 % ceil(par.theta_states/2)% theta goal reached
        if sP(2)==ceil(par.q_states/2) % q goal reached
            t = 1; %termiante
        else 
            t=0;   %continue
        end
    else
        t=0;       %continue
    end  
end

function a = execute_policy(Q, s, par)
    % Selects an action for state s using the
    % epsilon-greedy algorithm.
    
    %epsilon-greedy
    prob = rand;
    if prob <= par.epsilon %random action (explore)
        a = randi(par.actions);  
    else                   %current best action (exploit, greedy)
        %arg max_a Q(s,a)
        [~,I] = max(Q(s(1),s(2),:));
        a = I; 
    end
end

function Q = update_Q(Q, s, a, r, sP, aP, par)
    % SARSA update rule.
    Q(s(1),s(2),a) = Q(s(1),s(2),a)+ par.alpha*(r+ par.gamma*Q(sP(1),sP(2),aP)-Q(s(1),s(2),a));
end











