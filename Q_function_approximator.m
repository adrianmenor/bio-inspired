%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code approximates the learned Q function over the continuous
% state-space.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%   Clearing workspace
clear all
close all
clc

%%  Loading learned Q
load learned_Q Q;

% Plotting Q 
% subplot(2,1   ,1);
Q_plotter(Q,15);
% subplot(2,1,2);
% Q_plotter(Q,6);

% Plotting policy
% policy = policy_plotter(Q);
% save('learned_policy.mat','policy')

%%  Functions
function Z = Q_plotter(Q,u)
    XX = size(Q(:,:,u),1);
    YY = size(Q(:,:,u),2);
    [X,Y] = meshgrid(1:1:XX,1:1:YY);
    Z = Q(:,:,u)';
    surf(X,Y,Z);
    view(45,45)
    title('Q-values of discrete action \delta_e =',u,'FontSize',40);
    xlabel('\theta [-]','FontSize',40);
    ylabel('q [-]','FontSize',40);
    zlabel('Q-value','FontSize',40);
end

function policy = policy_plotter(Q)
    theta_size = size(Q(:,1,1),1);
    q_size = size(Q(1,:,1),2);
    policy = zeros(theta_size,q_size);
    
    for i = 1:theta_size % theta states
        for j =  1:q_size % q states
             %arg max_a Q(s,a)
            [~,I] = max(Q(i,j,:));
            policy(i,j) = I;
        end
    end
    [X,Y] = meshgrid(1:1:theta_size,1:1:q_size);
    Z = policy;
    surf(X,Y,Z);
    view(0,90)
    ax.FontSize = 100;
    c = colorbar
    c.Label.String = 'Discrete \delta_e input';
    
    
    title('Policy function \pi','FontSize',40 );
    xlabel('\theta [-]','FontSize',40);
    ylabel('q [-]','FontSize',40);
    zlabel('\delta_e','FontSize',40);
end


