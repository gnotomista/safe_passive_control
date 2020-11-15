% A Safety and Passivity Filter for Robot Teleoperation Systems
% Gennaro Notomista, Xiaoyi Cai
% International Workshop on Human-Friendly Robotics (HFR) 2020

clc
clear
close all

% system state (position and velocity)
% x = [p_x
%      p_y
%      v_x
%      v_y]
% system control input (acceleration)
% u = [a_x
%      a_y]
% system output (velocity)
% y = [v_x
%      v_y]

% flags
CASE = 'nothing'; % options: 'nothing', 'passive', 'safe', 'passafe' (=passive+safe)

% constants
DT = 0.002;
T = 10;

% system parameters, control gains, safe set parameters
SIGMA = 1;
K_P = 10;
K_D = 5;
K_I = 100;
D = 1;
GAMMA_1 = 1e1;
GAMMA_2 = 1e1;

% optimization parameters
QUADPROGOPTIONS = optimoptions(@quadprog,'Display','off');
U_MAX = inf;

% anonymous functions
gamma_u = @(s) (s>0)*1e2*s + (s<=0)*1e3*s;
gamma_x = @(s) 1e2*s^3;
robot_dynamics = @(x,u) [[zeros(2) eye(2);...
    zeros(2) -SIGMA*eye(2)]*x + [zeros(2);eye(2)]*u; ...
    x(3:4)];
uh = @(t,uh_prev) [-1; 0]*0.03 + uh_prev*0.99;
uh_dot = @(t,uh_prev) (uh(t,uh_prev)-uh_prev)/DT;
u_hat = @(x,t,uh_prev) -K_P*x(1:2)-K_D*x(3:4) + uh(t,uh_prev);
phi = @(x,u,t,uh_prev) -K_P*x(3:4)-K_D*(-SIGMA*x(3:4)+u)+K_I*(u_hat(x,t,uh_prev)-u)+uh_dot(t,uh_prev);
Au = @(x) x(3:4)';
bu = @(x,u,t,uh_prev) -(1+3*SIGMA^2)*norm(x(3:4))^2+2*SIGMA*x(1:2)'*x(3:4)-(2*x(1:2)'-3*SIGMA*x(3:4)')*u-x(3:4)'*phi(x,u,t,uh_dot(t,uh_prev))+gamma_u(2*SIGMA*norm(x(3:4))^2-2*x(1:2)'*x(3:4)-x(3:4)'*u);
Ax = @(x) -2*x(1:2)';
bx = @(x,u,t,uh_prev) (2*GAMMA_1*GAMMA_2*x(1:2)'+2*(GAMMA_1+GAMMA_2-SIGMA)*x(3:4)'+2*u')*x(3:4)+(4*x(3:4)'+2*(GAMMA_1+GAMMA_2-SIGMA)*x(1:2)')*(-SIGMA*x(3:4)+u)+2*x(1:2)'*phi(x,u,t,uh_prev)+gamma_x(GAMMA_1*GAMMA_2*norm(x(1:2))^2+2*norm(x(3:4))^2+2*(GAMMA_1+GAMMA_2-SIGMA)*x(1:2)'*x(3:4)+2*x(1:2)'*u-GAMMA_1*GAMMA_2*D^2);
V = @(x) norm(x)^2;

% initialize system variables
x = [2;2;-1;-1];
u = [0;0];
uh_prev = [0;0];
V_prev = V(x);

% initialize plot
figure('units','pixels','position',[0 0 1920 1080]), hold on, grid on, grid minor, axis equal, axis([-3 3 -3 3])
h_robot = scatter(x(1), x(2), 1000, '.');
h_robot_traj = line(x(1), x(2), 'LineStyle', '--', 'LineWidth', 2);
h_safe = plot(D*cos(linspace(0,2*pi,100)), D*sin(linspace(0,2*pi,100)), 'LineWidth', 2); uistack(h_safe, 'bottom')
h_safe_patch = patch(D*cos(linspace(0,2*pi,100)), D*sin(linspace(0,2*pi,100)), 'r', 'FaceAlpha', 0.1); uistack(h_safe_patch, 'bottom')
scatter(0, 0, 1000, '.');

% initialize log variables
robot_state_trajectory = nan(4,length(0:DT:T));
robot_input_trajectory = nan(2,length(0:DT:T));
power_margin = nan(1,length(0:DT:T));
safety_margin = nan(1,length(0:DT:T));

for t = 0 : DT : T
    n = round(t/DT+1);
    
    % evaluate controller
    if strcmp(CASE, 'nothing')
        Aqp = [];
        Bqp = [];
    elseif strcmp(CASE, 'passive')
        Aqp = Au(x);
        Bqp = bu(x,u,t,uh_prev);
    elseif strcmp(CASE, 'safe')
        Aqp = Ax(x);
        Bqp = bx(x,u,t,uh_prev);
    elseif strcmp(CASE, 'passafe')
        Aqp = [Au(x); Ax(x)];
        Bqp = [bu(x,u,t,uh_prev); bx(x,u,t,uh_prev)];
    end
    v_star = quadprog(eye(2), zeros(2,1), Aqp, Bqp, [], [], -U_MAX*ones(2,1), U_MAX*ones(2,1), [], QUADPROGOPTIONS);
    
    % dynamic equations
    x_dot_y = robot_dynamics(x,u);
    x_dot = x_dot_y(1:4);
    y = x_dot_y(5:6);
    u_dot = phi(x,u,t,uh_prev) + v_star;
    
    % log state and input trajectories, and power and safety margin
    robot_state_trajectory(:,n) = x;
    robot_input_trajectory(:,n) = u;
    power_margin(n) = u'*y - 2*x'*x_dot;
    safety_margin(n) = norm(x(1:2))^2-1;
    
    % dynamics simulation step
    x = x + x_dot * DT;
    u = u + u_dot * DT;
    
    % update previous values of V and uh
    V_prev = V(x);
    uh_prev = uh(t,uh_prev);
    
    % update plots
    h_robot.XData = x(1);
    h_robot.YData = x(2);
    h_robot_traj.XData = robot_state_trajectory(1,1:10:end);
    h_robot_traj.YData = robot_state_trajectory(2,1:10:end);
    drawnow limitrate
end

% passivity and safety margin plots
figure, hold on, grid on, grid minor, axis tight
plot(power_margin, 'LineWidth', 2);
plot(safety_margin, 'LineWidth', 2);
legend('Power margin', 'Safety margin')