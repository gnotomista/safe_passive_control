% Passivity-Based Decentralized Control of Multi-Robot Systems With Delays Using Control Barrier Functions
% Gennaro Notomista, Xiaoyi Cai, Junya Yamauchi, Magnus Egerstedt
% International Symposium on Multi-Robot and Multi-Agent Systems (MRS) 2019

clc
clear
close all

addpath('include')
run('swarm_sim/init.m')

% flags
USE_CBFS = true;

% anonymous functions
alpha = @(x) x; % class K function

% simulation parameters
T_MAX = 1e3;
DT = 0.033; % simulation time step

% graph Laplacian
L = [3 -1 -1 -1 0 -1 ; ...
    -1 3 -1 0 -1 0 ; ...
    -1 -1 3 -1 0 -1 ; ...
    -1 0 -1 3 -1 0 ; ...
    0 -1 0 -1 3 -1 ; ...
    -1 0 -1 0 -1 3];
N = size(L,1); % number of robots
% weight matrix
D = 0.4;
D_DIAG = 2*D;
WEIGHTS = [0 D sqrt(3)*D D_DIAG 0 D; ...
    D 0 D 0 D_DIAG 0; ...
    sqrt(3)*D D 0 D 0 D_DIAG; ...
    D_DIAG 0 D 0 D 0; ...
    0 D_DIAG 0 D 0 D; ...
    D 0 D_DIAG 0 D 0];

% PD control gains
K_P = 5; % proportional
K_D = 1; % derivative

% optimization parameters
QUADPROGOPTIONS = optimoptions(@quadprog,'Display','off');
RHO = 1e-9;
U_MAX = inf;

% plot settings
AXIS_LIM = 2;

% multi-robot system (initialize N robots modeled with unicycle dynamics)
ROBOT_W = 0.1;
ROBOT_L = 0.1;
robots = cell(1,N);
for i = 1 : N
    for i = 1 : N
        robots{i} = Unicycle('width', ROBOT_W, ...
            'length', ROBOT_L, ...
            'initialState', [-AXIS_LIM/2+AXIS_LIM*rand(2,1);2*pi*rand()], ...
            'simulationTimeStep', DT);
    end
end
s = Swarm('robots', robots, 'L', L);
x = s.getPoses();

% initialize fixed-length queues to simulate delays
q = cell(N,N);
for i = 1 : N
    for j = i+1 : N
        delay = randi(11)-1; % random delay between 0 and 10 iterations
        q{i,j} = FixedLengthQueue(1+delay);
        for n = 1 : 1+delay
            q{i,j}.add([x(1:2,j);zeros(2,1)]);
        end
    end
end
for j = 1 : N
    for i = j+1 : N
        % delay = q{j,i}.l-1; % symmetric delays
        delay = randi(11)-1; % non-symmetric delays: different delays in exchanging position information for i->j and j->i
        q{i,j} = FixedLengthQueue(1+delay);
        for n = 1 : 1+delay
            q{i,j}.add([x(1:2,j);zeros(2,1)]);
        end
    end
end
for i = 1 : N
    q{i,i} = FixedLengthQueue(1);
    q{i,i}.add([x(1:2,i);zeros(2,1)]);
end

% initialize CBFs values
h = zeros(1,N);

% figure
s.plotFigure()
axis(AXIS_LIM*[-1 1 -1 1])

% main loop
for t = 1 : T_MAX
    x = s.getPoses();
    
    % evaluate "passive" formation controller
    dx = zeros(2,N);
    for i = 1 : N
        % evaluate nominal formation controller (PD for damped double integrator dynamics)
        ui_hat = zeros(2,1);
        xi = q{i,i}.peek();
        xi(1:2) = x(1:2,i);
        for j = s.getNeighbors(i)
            xj = q{i,j}.peek();
            ui_hat = ui_hat + K_P * (norm(xj(1:2)-xi(1:2))^2-WEIGHTS(i,j)^2) * (xj(1:2)-xi(1:2));
        end
        % if USE_CBFS, modify the nominal formation controller to ensure a passivity condition in integral form
        if USE_CBFS
            A = [xi(3:4)' -norm(xi(3:4))^2];
            b = alpha(h(i))+K_D*norm(xi(3:4))^2;
            ui_a_star = quadprog(blkdiag(2*eye(2),2*RHO), [-2*ui_hat' 0], A, b, [], [], [-U_MAX -U_MAX 0], [U_MAX U_MAX inf], [], QUADPROGOPTIONS);
            ui_star = ui_a_star(1:2);
            Kd_star = K_D + ui_a_star(3);
        else
            ui_a_star = zeros(3,1);
            ui_star = min(U_MAX, max(-U_MAX, ui_hat));
            Kd_star = K_D;
        end
        % update CBFs values
        h(i) = h(i) + (-ui_hat'*xi(3:4)+Kd_star*norm(xi(3:4))^2) * DT;
        % integrate the damped double integrator dynamics to calculate the velocity input dx for robot i
        xi = xi + (kron([0 1;0 -Kd_star],eye(2))*xi+kron([0;1],eye(2))*ui_star) * DT;
        dx(:,i) = xi(3:4);
    end
    
    % enqueue robot poses into communication channels
    for i = 1 : N
        q{i,i}.add([x(1:2,i);dx(:,i)]);
        for j = s.getNeighbors(i)
            q{j,i}.add([x(1:2,i);dx(:,i)]);
        end
    end
    
    % dynamics simulation step
    s.moveSingleIntegrators(dx);
    
    % plot robots and graph corresponding to the Laplacian L
    s.plotRobots()
    s.plotGraph()
    drawnow limitrate
end
