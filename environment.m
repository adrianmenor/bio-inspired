%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code simulates the rocket/aircraft dynamics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function xe = environment(x,u,par)
%     % System ss model
% %     A = [1 1; 1 0];
%      A = [0 1; 14.7805 0];
%     B = [0; 3.4858];
%    
%     % Generating new state
%     x_dot = A*x'+B*u;
%     xe = (x'+x_dot*par.simstep)';
% end

% plane-pitch environment
function xe = environment(x,u,par)
    % System ss model
    A = [-2.114 0; 1 0];
    B = [-12.64; 0];
   
    % Generating new state
    x_dot = A*x'+B*u;
    xe = (x'+x_dot*par.simstep)';
end