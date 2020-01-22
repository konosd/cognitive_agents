% parameters, number of agents, trajectories, etc.
n_agent = 200;       %number of agents
n_vsteps = 20;      %number of virtual steps
n_steps = 1000;       %number of real steps
n_traj = 20;        %number of trajectories
sigma = 5;          %diameter
box_length = 80*sigma;    %area explored
timestep = 0.1;     % dt timestep
friction = 1.0;     %gamma
temperature = 20.0;  %temperature, to test if ss achieved
noise = sqrt(2.0*friction*temperature/timestep);
repul_strength = 20.0;
repul_exp = 60.0;



% ------------- Initialization-- ------------------------------------------
agent_velo = zeros(n_agent,2);
agent_coor = initialize_agents(n_agent, sigma, box_length);
force_init = hard_repulsion(agent_coor, sigma, box_length, repul_strength);
%plot_agents(agent_coor, sigma, box_length, force_init, 1)

% -------------------------------------------------------------------------

% Initialize coordinates
everything_coor_x = zeros(n_agent, n_steps);
everything_coor_y = zeros(n_agent, n_steps);
everything_coor_x(:,1) = agent_coor(:,1);
everything_coor_y(:,1) = agent_coor(:,2);

% first timestep
all_gyrations = calc_all_gyrations(n_agent, n_traj, agent_coor,...
    agent_velo, sigma, box_length, repul_strength, friction, noise,...
    timestep, n_vsteps);

[new_coor, new_velo] = next_timestep( agent_coor, ...
    agent_coor, agent_velo, timestep, ...
    n_agent, all_gyrations, n_traj, sigma, box_length, repul_strength, ...
    friction, noise);
everything_coor_x(:,2) = mod(new_coor(:,1),box_length);
everything_coor_y(:,2) = mod(new_coor(:,2), box_length);

%next timesteps
for step=3:n_steps
    agent_coor = [everything_coor_x(:,step-1) everything_coor_y(:,step-1)];
    past_agent_coor = [everything_coor_x(:,step-2) everything_coor_y(:,step-2)];
    all_gyrations = calc_all_gyrations(n_agent, n_traj, agent_coor,...
        agent_velo, sigma, box_length, repul_strength, friction, noise,...
        timestep, n_vsteps);

    [new_coor, new_velo] = next_timestep( agent_coor, ...
        past_agent_coor, agent_velo, timestep, ...
        n_agent, all_gyrations, n_traj, sigma, box_length, repul_strength, ...
        friction, noise);
    everything_coor_x(:,step) = mod(new_coor(:,1), box_length);
    everything_coor_y(:,step) = mod(new_coor(:,2), box_length);
    %plot_agents(agent_coor, sigma, box_length, force_init, 1)
end

% plot the real trajectories:
% for agent = 1:n_agent
%     traj = [everything_coor_x(agent,:); everything_coor_y(agent,:)].';
%     plot_trajectory(traj, box_length, rand(1,3))
% end

if true
    vobj = VideoWriter('simulation.avi', 'MPEG-4');
    vobj.FrameRate = 10;
    vobj.open;
    fig = figure('visible','off');
    for step = 1:n_steps
        plot_coor = [everything_coor_x(:,step) everything_coor_y(:,step)];
        force_repulsive = hard_repulsion(plot_coor, sigma, box_length, ...
            repul_strength);
        plot_agents(plot_coor, sigma, box_length, force_repulsive, 1);
        frame = getframe(fig);
        writeVideo(vobj,frame);
    end
    vobj.close;
    vobj.delete;
end


% ---------------------- End of main --------------------------------------
% -------------------------------------------------------------------------
% --------------------- Helper Functions ----------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------- Sample from Bivariate Distribution -----------------
% -------------------------------------------------------------------------
function bivar_random = bivariate_normal( samples )
    bivar_random = repmat( [0 0], samples, 1) + ...
        randn(samples, 2) * [1 0; 0 1];
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------- Virtual Trajectory Creation ------------------------
% -------------------------------------------------------------------------
% i is the agent we are investigating

function [traj_coor, traj_velo] = virtual_traj(agent_coor, agent_velo, i, ...
    n_virt_steps, diameter, area, strength, friction, noise, dt)
    
    % First element is the real one, so updating number of virtual time
    % steps
    n_virt_steps = n_virt_steps + 1;
    
    % Initialize vectors of path and velocities
    traj_coor = zeros(n_virt_steps,2);
    traj_velo = zeros(n_virt_steps,2);
    grid_coor = agent_coor;
    
    % Generate random numbers for noise for every step
    virt_steps_noise = noise * bivariate_normal(n_virt_steps);
    
    
    % starting the iteration for the first virtual timestep
    traj_coor(1,:) = agent_coor(i,:);
    traj_velo(1,:) = agent_velo(i,:);
    
    f_rep = hard_repulsion_agent(agent_coor, i, diameter, area, strength);
    f_langevin = -friction * traj_velo(1,:) + virt_steps_noise(1,:);
    f_tot = f_langevin + f_rep;
    traj_coor(2,:) = traj_coor(1,:) + traj_velo(1,:)* dt + ...
        0.5 * f_tot * (dt^2) ;
    traj_velo(2,:) = traj_velo(1,:) + f_tot* dt ;
    
    % Update the grid;
    grid_coor(i,:) = traj_coor(2,:);
    
    for j=3:n_virt_steps
        % find repulsion force for step
        f_rep = hard_repulsion_agent(grid_coor, i, diameter, area, strength);
        f_langevin = -friction * traj_velo(j-1,:) + virt_steps_noise(j,:);
        f_tot = f_rep + f_langevin;
        % update velocity and position of virtual timestep
        traj_velo(j,:) = traj_velo(j-1,:) + f_tot * dt;
        traj_coor(j,:) = 2 * traj_coor(j-1,:) - traj_coor(j-2,:) + ...
            f_tot * (dt^2);
        % update the grid:
        traj_coor(j,:) = mod(traj_coor(j,:),area);
        grid_coor(i,:) = traj_coor(j,:);
    end
    
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% ------------- Find gyration for  trajectory------------------------------
% -------------------------------------------------------------------------
function gyration = calc_gyration(trajectory)
    mean_gyration = mean(trajectory, 1);
    gyration = zeros(1,2);
    for i=2:length(trajectory)
        gyration = gyration + (trajectory(i,:) - mean_gyration).^2 ;
    end
    gyration = gyration/(length(trajectory)-1);
    gyration = sqrt(gyration(1)^2 + gyration(2)^2);

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
            Dx = agent_coordinates(i,1) - agent_coordinates(j,1);
            Dy = agent_coordinates(i,2) - agent_coordinates(j,2);
            if Dx > 0.5 * area
                Dx = -(area + agent_coordinates(j,1) - agent_coordinates(i,1));
            elseif Dx < -0.5 *area
                Dx = area - agent_coordinates(j,1) + agent_coordinates(i,1);
            end
            if Dy > 0.5 * area
                Dy = -(area + agent_coordinates(j,2) - agent_coordinates(i,2));
            elseif Dy < -0.5 *area
                Dy = area - agent_coordinates(j,2) + agent_coordinates(i,2);
            end
            ag_dist = sqrt(Dx^2 + Dy^2);
            if i ~= j && ag_dist < 2 * diameter
                magnitude = strength*(2*diameter - ag_dist);
                force_rep(i,1) = force_rep(i,1) + magnitude * Dx/ ag_dist;
                force_rep(i,2) = force_rep(i,2) + magnitude * Dy/ ag_dist;
            end
        end
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% ------------- Hardcore Repulsion for One Agent---------------------------
% -------------------------------------------------------------------------
function force_rep = hard_repulsion_agent(agent_coordinates, i, ...
        diameter, area, strength)
    force_rep = zeros(1, 2);
    for j = 1:length(agent_coordinates)
        Dx = agent_coordinates(i,1) - agent_coordinates(j,1);
        Dy = agent_coordinates(i,2) - agent_coordinates(j,2);
        if Dx > 0.5 * area
            Dx = -(area + agent_coordinates(j,1) - agent_coordinates(i,1));
        elseif Dx < -0.5 *area
            Dx = area - agent_coordinates(j,1) + agent_coordinates(i,1);
        end
        if Dy > 0.5 * area
            Dy = -(area + agent_coordinates(j,2) - agent_coordinates(i,2));
        elseif Dy < -0.5 *area
            Dy = area - agent_coordinates(j,2) + agent_coordinates(i,2);
        end
        ag_dist = sqrt(Dx^2 + Dy^2);
        if i ~= j && ag_dist < 2 * diameter
            magnitude = strength*(2*diameter - ag_dist);
            force_rep(1) = force_rep(1) + magnitude * Dx/ ag_dist;
            force_rep(2) = force_rep(2) + magnitude * Dy/ ag_dist;
        end
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
               

% -------------------------------------------------------------------------
% -------------------- Calculate All Gyrations ----------------------------
% -------------------------------------------------------------------------
function all_gyrations = calc_all_gyrations(n_agent, n_traj, agent_coor,...
    agent_velo, sigma, box_length, repul_strength, friction, noise,...
    timestep, n_vsteps)

    all_gyrations = zeros(n_agent, n_traj);

    for agent_no = 1:n_agent
        %disp('Entering agent')
        for tr=1:n_traj

            [virt_coor, virt_velo] = virtual_traj(agent_coor, agent_velo,...
                agent_no, n_vsteps, ...
            sigma, box_length, repul_strength, friction, noise, timestep);

            %plot_trajectory(virt_coor,box_length, rand(1,3))
            all_gyrations(agent_no, tr) = calc_gyration(virt_coor);
        end
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------- One Timestep calculation ---------------------------
% -------------------------------------------------------------------------
function [new_coor, new_velo] = next_timestep( agent_coor, ...
    past_agent_coor, agent_velo, timestep, ...
    n_agent, all_gyrations, n_traj, diameter, area, strength, friction, noise)
    new_coor = zeros(n_agent, 2);
    new_velo = zeros(n_agent, 2);
    for agent=1:n_agent
        % Calculate langevin force based on gyration
        f_lang_agent = zeros(1,2);
        agent_mean_gyration = mean(all_gyrations(agent,:));
        for j=1:n_traj
            norm_gyr = log(all_gyrations(agent,j)/agent_mean_gyration);
            f_lang_agent = f_lang_agent + norm_gyr * bivariate_normal(1);
        end
        f_lang_agent = f_lang_agent/n_traj;
        f_lang_agent = -friction * agent_velo(agent,:) + ...
            noise * f_lang_agent;
        f_rep = hard_repulsion_agent(agent_coor, agent, diameter, ...
            area, strength);
        f_tot = f_lang_agent + f_rep;
        new_coor(agent, :) = 2 * agent_coor(agent, :) - ...
            past_agent_coor(agent,:) + f_tot * (timestep^2);
        new_velo(agent, :) = agent_velo(agent, :) + timestep * f_tot;
        
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------





% -------------------------------------------------------------------------
% ------------- Function to plot screenshot--------------------------------
% -------------------------------------------------------------------------
function plot_agents(agent_coordinates, diameter, box_length, forces, scl)
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
    quiver(agent_coordinates(:,1), agent_coordinates(:,2), ...
        forces(:,1), forces(:,2), scl);
    grid on;
    xlim([0, box_length]);
    ylim([0, box_length])
    hold off;
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% --------------------Plot Virtual Trajectory -----------------------------
% -------------------------------------------------------------------------
function plot_trajectory(traj_coord, box_length, rand_color)
    counter = 1;
    curve = animatedline('linewidth',1, 'color', rand_color);
    for i=1:(length(traj_coord)-1)
        dif = norm(traj_coord(i+1,:)-traj_coord(i,:));
        if dif > box_length * 0.5
            addpoints(curve, traj_coord(counter:i,1), traj_coord(counter:i,2));
            drawnow
            curve = animatedline('linewidth',1, 'color', rand_color);
            counter = i+1;
        end
    end
    addpoints(curve, traj_coord(counter:end,1), traj_coord(counter:end,2));
    drawnow   
end