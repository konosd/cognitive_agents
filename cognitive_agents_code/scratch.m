% parameters, number of agents, trajectories, etc.
n_agent = 20;       %number of agents
n_vsteps = 20;      %number of virtual steps
n_steps = 20;       %number of real steps
n_traj = 20;        %number of trajectories
sigma = 5;          %diameter
box_length = 80;    %area explored
timestep = 0.1;     % dt timestep
friction = 1.0;     %gamma
temperature = 1.0;  %temperature, to test if ss achieved
noise = sqrt(2.0*friction*temperature/timestep);
repul_strength = 20.0;
repul_exp = 60.0;


agent_velo = zeros(n_agent,2);
agent_coor = initialize_agents(n_agent, sigma, box_length);

force_init = hard_repulsion(agent_coor, sigma, box_length, repul_strength);


plot_agents(agent_coor, sigma, box_length)

% ---------------------- End of main --------------------------------------
% -------------------------------------------------------------------------
% ------------- Helper Functions ------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% ------------- Function to plot screenshot--------------------------------
% -------------------------------------------------------------------------
function plot_agents(agent_coordinates, diameter, box_length)
    theta = 0:0.01:2*pi;
    for k = 1: length(agent_coordinates)
        xCenter = agent_coordinates(k,1);
        yCenter = agent_coordinates(k,2);
        thisX = diameter * 0.5 * cos(theta) + xCenter;
        thisY = diameter * 0.5 * sin(theta) + yCenter;
        plot(xCenter, yCenter, 'r+', 'MarkerSize', 10, 'LineWidth', 1);
        text(xCenter, yCenter, num2str(k));
        hold on;
        plot(thisX, thisY, 'b-', 'LineWidth', 2);
    end
    grid on;
    xlim([0, box_length]);
    ylim([0,box_length])
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% ------------- Function to initialize agents------------------------------
% -------------------------------------------------------------------------
function init_agents = initialize_agents(n_agents, agent_diameter, box_dim)
    
    init_agents = zeros(n_agents, 2);

    % --------- First set randomly the agents in the box --------------
    % make sure we have enough space
    theor_max = (fix(box_dim*0.5))^2;     %maximum number of agents
    min_distance = 1.1 * agent_diameter;
    radius = agent_diameter * 0.5;
    allowed_length = box_dim - agent_diameter;
    if n_agents >= theor_max
        disp("theoretical maximum limit exceeded ");
    end

    if n_agents >= theor_max * 2/3;
        disp("Slow regime limit exceeded ");
    end

    init_agents(1,:) = allowed_length *rand(2,1) + radius;

    % then position everyone else:

    for i = 2:n_agents
        next_found = 0;
        while next_found == 0
            candidate_coor = allowed_length * rand(2,1) + radius;
            all_far = 0;
            for j = 1:i-1
                dist_vector = norm(init_agents(j,:) - candidate_coor.');
                all_far = all_far + (dist_vector < min_distance);
            end
            if all_far == 0
                init_agents(i,:) = candidate_coor;
                next_found = 1;
            end
        end
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% ------------- Hardcore Repulsion ----------------------------------------
% -------------------------------------------------------------------------
function force_rep = hard_repulsion(agent_coordinates, diameter, area, strength)
    force_rep = zeros(length(agent_coordinates), 2);
    for i = 1:length(agent_coordinates)
        for j = 1:length(agent_coordinates)
            x_dist = norm(agent_coordinates(i,1) - agent_coordinates(j,1));
            y_dist = norm(agent_coordinates(i,2) - agent_coordinates(j,2));
            if x_dist > 0.5 * area
                x_dist = area - x_dist;
            end
            if y_dist > 0.5 * area
                y_dist = area - y_dist;
            end
            ag_dist = sqrt(x_dist^2 + y_dist^2);
            if i ~= j && ag_dist < 2 * diameter
                theta = atan((agent_coordinates(i,2)-agent_coordinates(j,2))/ ...
                    (agent_coordinates(i,1)-agent_coordinates(j,1)));
                magnitude = strength*(2*diameter - ag_dist)
                force_rep(i,1) = force_rep(i,1) + magnitude * ...
                    (agent_coordinates(i,1)-agent_coordinates(j,1))/ ...
                    ag_dist;
                force_rep(i,2) = force_rep(i,2) + magnitude * ...
                    (agent_coordinates(i,2)-agent_coordinates(j,2))/ ...
                    ag_dist;
            end
        end
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
               