%%
% parameters, number of agents, trajectories, etc.
n_agent = 2;       %number of agents
n_vsteps = 20000;      %number of virtual steps
n_steps = 1;       %number of real steps
n_traj = 20;        %number of trajectories
sigma = 1;          %diameter
box_length = 80*sigma;    %area explored
timestep = 0.1;     % dt timestep
friction = 0.2;     %gamma
temperature = 1.0;  %temperature
noise = sqrt(2.0*friction*temperature/timestep);
repul_strength = 20.0;
repul_exp = 60.0;
repul_type = "soft";
v_repul_type = "hard";
pi = 4 * atan(1);

% Add synthetic agents - this is where you define whether an agent will be
% able to move or not. If the number of the agent is in the synthetic
% vector then the agent will not move
synthetic = [];


% parameters for the active brownian agens. Additional ones are: gamma(r)
% additive friction, U(r) potential field, or modify q(r) as an intake
% field. Here, q is only constant. Noise is the same as with passive BP.
q0 = 1.0;    % energy intake from the environment
food_radius = 10;
food_center = [80*sigma*0.5 80*sigma*0.5];
q = @(x1, x2) q0 * ( ((x1-food_center(1))^2 + (x2-food_center(2))^2) < food_radius^2 );
d2 = 1.0;   % conversion rate of internal-to-kinetic energy
c = 0.1;    % dissipation of internal energy
abp_steps = 60000;
h = 0.1;
D = 0.001;

% Potential field
a = 0;
U = @(x1,x2) 0.5 * a * (x1^2 + x2^2);


% Filling fraction:
phi = n_agent * pi * sigma^2 / (4* box_length^2);
disp("Filling fraction is " + phi)

%% ------------- Initialization--------------------------------------------
agent_coor = initialize_agents(n_agent, sigma, box_length);
agent_velo = zeros(n_agent,2);
% For the active particles, add energy depot
agent_energies = zeros(n_agent,1);

% % To start from a previous step
% agent_coor = [all_x(:,1000) all_y(: , 1000)];
% agent_velo = [vel_x(:,1000) vel_y(: , 1000)];

force_init = repulsion(agent_coor, sigma, box_length, repul_strength, repul_type);


% -------------------------------------------------------------------------



%% ---------------- Recreating model ABP  paper----------------------------

% % Initial positions
% e0 = 0;
% r0 = [1 0];
% v0 = [1 0.5];
% 
% % Define deterministic ODEs of motion of an active brownian agent
% f = {@(v, e) v;
%     @(v, e) (d2 * e *v - friction * v);
%     @(v, e) q - c * e - d2 * (v(1)^2 + v(2)^2) * e};
% 
% [r,v,e] = rk4(f, r0, v0, e0, d2, friction, q, c, h, abp_steps);
% 
% t = [0:h:abp_steps*h];
% 
% % To reproduce figure 3 of abp_model paper:
% % Initial positions
% x1 = zeros(abp_steps+1,1);
% x2 = zeros(abp_steps+1,1);
% v1 = zeros(abp_steps+1, 1);
% v2 = zeros(abp_steps+1, 1);
% e = zeros(abp_steps+1, 1);
% 
% x1(1) = 0; 
% x2(1) = 0;
% v1(1) = 0;
% v2(1) = 0;
% 
% f = {@(x1, x2, v1, v2, e) v1;
%     @(x1, x2, v1, v2, e) v2;
%     @(x1, x2, v1, v2, e) d2 * e *v1 - friction * v1 - a * x1 ;
%     @(x1, x2, v1, v2, e) d2 * e *v2 - friction * v2 - a * x2;
%     @(x1, x2, v1, v2, e) q(x1, x2) - c* e - d2* e* (v1^2 + v2^2) };
%     
% % RK4 integration
% for i=1:abp_steps
% 
%     k1 = h * [ f{1}(x1(i), x2(i), v1(i), v2(i), e(i));
%         f{2}(x1(i), x2(i), v1(i), v2(i), e(i));
%         f{3}(x1(i), x2(i), v1(i), v2(i), e(i));
%         f{4}(x1(i), x2(i), v1(i), v2(i), e(i));
%         f{5}(x1(i), x2(i), v1(i), v2(i), e(i))];
%     k2 = h * [ f{1}(x1(i) + k1(1)/2, x2(i)+ k1(2)/2, v1(i)+ k1(3)/2, v2(i)+ k1(4)/2, e(i)+ k1(5)/2);
%         f{2}(x1(i) + k1(1)/2, x2(i)+ k1(2)/2, v1(i)+ k1(3)/2, v2(i)+ k1(4)/2, e(i)+ k1(5)/2);
%         f{3}(x1(i) + k1(1)/2, x2(i)+ k1(2)/2, v1(i)+ k1(3)/2, v2(i)+ k1(4)/2, e(i)+ k1(5)/2);
%         f{4}(x1(i) + k1(1)/2, x2(i)+ k1(2)/2, v1(i)+ k1(3)/2, v2(i)+ k1(4)/2, e(i)+ k1(5)/2);
%         f{5}(x1(i) + k1(1)/2, x2(i)+ k1(2)/2, v1(i)+ k1(3)/2, v2(i)+ k1(4)/2, e(i)+ k1(5)/2)];
%     k3 = h * [ f{1}(x1(i) + k2(1)/2, x2(i)+ k2(2)/2, v1(i)+ k2(3)/2, v2(i)+ k2(4)/2, e(i)+ k2(5)/2);
%         f{2}(x1(i) + k2(1)/2, x2(i)+ k2(2)/2, v1(i)+ k2(3)/2, v2(i)+ k2(4)/2, e(i)+ k2(5)/2);
%         f{3}(x1(i) + k2(1)/2, x2(i)+ k2(2)/2, v1(i)+ k2(3)/2, v2(i)+ k2(4)/2, e(i)+ k2(5)/2);
%         f{4}(x1(i) + k2(1)/2, x2(i)+ k2(2)/2, v1(i)+ k2(3)/2, v2(i)+ k2(4)/2, e(i)+ k2(5)/2);
%         f{5}(x1(i) + k2(1)/2, x2(i)+ k2(2)/2, v1(i)+ k2(3)/2, v2(i)+ k2(4)/2, e(i)+ k2(5)/2)];
%     k4 = h * [ f{1}(x1(i) + k3(1), x2(i)+ k3(2), v1(i)+ k3(3), v2(i)+ k3(4), e(i)+ k3(5));
%         f{2}(x1(i) + k3(1), x2(i)+ k3(2), v1(i)+ k3(3), v2(i)+ k3(4), e(i)+ k3(5));
%         f{3}(x1(i) + k3(1), x2(i)+ k3(2), v1(i)+ k3(3), v2(i)+ k3(4), e(i)+ k3(5));
%         f{4}(x1(i) + k3(1), x2(i)+ k3(2), v1(i)+ k3(3), v2(i)+ k3(4), e(i)+ k3(5));
%         f{5}(x1(i) + k3(1), x2(i)+ k3(2), v1(i)+ k3(3), v2(i)+ k3(4), e(i)+ k3(5))];
%     noise = sqrt(2*D)*bivariate_normal(1);
%     x1(i+1) = x1(i) + (k1(1) + 2*k2(1) + 2*k3(1) + k4(1))/6 ;
%     x2(i+1) = x2(i) + (k1(2) + 2*k2(2) + 2*k3(2) + k4(2))/6 ;
%     v1(i+1) = v1(i) + (k1(3) + 2*k2(3) + 2*k3(3) + k4(3))/6 + noise(1) ;
%     v2(i+1) = v2(i) + (k1(4) + 2*k2(4) + 2*k3(4) + k4(4))/6 + noise(2) ;
%     e(i+1) = e(i) + (k1(5) + 2*k2(5) + 2*k3(5) + k4(5))/6 ;
%     
% end
% 
% plot(x1, x2)


%% ---------- Collective ABP ----------------------------------------------
% what happens if i try to solve the problem with the same cognitive force,
% etc., but updated motions.

[my_traj_coor, my_traj_velo, bound, traj_init] = virtual_traj_abp(agent_coor, agent_velo, 1, ...
    n_vsteps, sigma, box_length, repul_strength, friction, D, h, v_repul_type, 0, d2, a, c, q);




fig = figure(1);
plot_agents(agent_coor, sigma, box_length, force_init, 1)
plot_trajectory(bound, box_length, rand(1,3))



%% ----------- Pipeline for passive BP, from old code  --------------------

% [my_traj_coor, my_traj_velo, bound, traj_init] = virtual_traj(agent_coor, agent_velo, 1, ...
%     n_vsteps, sigma, box_length, repul_strength, friction, noise, timestep, v_repul_type);
% 
% fig = figure(1);
% plot_agents(agent_coor, sigma, box_length, force_init, 1)
% plot_trajectory(bound, box_length, rand(1,3))

% ----------- To visualize all virtual trajectories of one agent ---------
% fig = figure(1);
% plot_agents(agent_coor, sigma, box_length, force_init, 1)
% for tr = 1:n_traj
%    [my_virt_traj, my_virt_velo, bound, traj_init] =  virtual_traj(agent_coor, agent_velo, 1, ...
%      n_vsteps, sigma, box_length, repul_strength, friction, noise, timestep, v_repul_type);
%     plot_trajectory(bound, box_length, rand(1,3))
% end

% ---------- To visualize all virtual trajectories of all agents ---------
% fig = figure(1);
% plot_agents(agent_coor, sigma, box_length, force_init, 1)
% for agent = 1:n_agent
%     for tr = 1:n_traj
%        [my_virt_traj, my_virt_velo, bound, traj_init] =  virtual_traj(agent_coor, ...
%             agent_velo, agent, ...
%          n_vsteps, sigma, box_length, repul_strength, friction, noise,...
%               timestep, v_repul_type);
%         plot_trajectory(bound, box_length, rand(1,3))
%     end
% end

% %% ---- To visualize everything, solve the full problem, chunk below ------
% % 
% [all_x, all_y, vel_x, vel_y] = cef_solver( agent_coor,...
%     n_agent, n_traj, sigma, ...
%     box_length, repul_strength, friction, noise, timestep, n_vsteps, ...
%     n_steps, repul_type, v_repul_type, false, synthetic);
% 
% incr = 4;
% coordat = zeros(n_agent * n_steps / incr, 4);
% for i=1:n_steps/4
%     coordat(((i-1)*n_agent+1):(i*n_agent) , :) = [all_x(:,i) all_y(:,i) ...
%         vel_x(:,i) vel_y(:,i)];
% end
% 
% save coor.dat coordat -ascii



%% ---------------------- End of main --------------------------------------

% -------------------------------------------------------------------------
% --------------------- Helper Functions ---------------------------------
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% --------------- Virtual traj for ABP   ----------------------------------
% -------------------------------------------------------------------------
function [traj_coor, traj_velo, ...
    bound_coor, traj_init_force] = virtual_traj_abp(agent_coor, agent_velo, i, ...
    n_virt_steps, diameter, area, strength, friction, D ,h, repul_type, e_init, d2, a, c, q)
    
    % First element is the real one, so updating number of virtual time
    % steps
    n_virt_steps = n_virt_steps;
    
    % Initialize vectors of path and velocities
    traj_coor = zeros(n_virt_steps+1,2);
    bound_coor = zeros(n_virt_steps+1,2);
    traj_velo = zeros(n_virt_steps+1,2);
    traj_energy = zeros(n_virt_steps+1,1);
    grid_coor = agent_coor;
    
    % Generate random numbers for noise for every step
    virt_steps_noise = bivariate_normal(n_virt_steps);
    % Save the first force to compute gyration later
    traj_init_force = virt_steps_noise(1,:);
    virt_steps_noise = sqrt(2*D) * virt_steps_noise;
    
    % The equations of motion and energy for the agent
    f = {@(x1, x2, v1, v2, e) v1;
    @(x1, x2, v1, v2, e) v2;
    @(x1, x2, v1, v2, e, fx) d2 * e *v1 - friction * v1 - a * x1 + fx;
    @(x1, x2, v1, v2, e, fy) d2 * e *v2 - friction * v2 - a * x2 + fy;
    @(x1, x2, v1, v2, e) q(x1, x2) - c* e - d2* e* (v1^2 + v2^2) };


    % starting the iteration for the first virtual timestep
    traj_coor(1,:) = agent_coor(i,:);
    traj_velo(1,:) = agent_velo(i,:);
    bound_coor(1,:) = agent_coor(i,:);
    traj_energy(1) = 0;
    
    
    for j=2:n_virt_steps
        [traj_coor(j,1), traj_coor(j,2), traj_velo(j,1),...
            traj_velo(j,2), traj_energy(j) ]= rk4_step(f, ...
            traj_coor(j-1,1), traj_coor(j-1,2), traj_velo(j-1,1),...
            traj_velo(j-1,2), traj_energy(j-1), h, diameter, ...
    area, strength, repul_type, virt_steps_noise(j,:), i, grid_coor);

        bound_coor(j,:) = mod(traj_coor(j,:), area);
        grid_coor(i,:) = mod(traj_coor(j,:), area);
    end
    
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% --------------------------- All timesteps -------------------------------
% -------------------------------------------------------------------------
function [everything_coor_x, everything_coor_y, all_velo_x, all_velo_y] =...
    cef_solver( agent_coor,...
    n_agent, n_traj, sigma, ...
    box_length, repul_strength, friction, noise, timestep, n_vsteps, ...
    n_steps, repul_type, v_repul_type, record, synthetic)
    % ------------- Initialization-- ------------------------------------------
    agent_velo = zeros(n_agent,2);
    force_init = repulsion(agent_coor, sigma, box_length, repul_strength, repul_type);
    %plot_agents(agent_coor, sigma, box_length, force_init, 1)

    % -------------------------------------------------------------------------

    % Initialize coordinates
    everything_coor_x = zeros(n_agent, n_steps);
    everything_coor_y = zeros(n_agent, n_steps);
    everything_coor_x(:,1) = agent_coor(:,1);
    everything_coor_y(:,1) = agent_coor(:,2);
    
    all_velo_x = zeros(n_agent, n_steps);
    all_velo_y = zeros(n_agent, n_steps);
    all_velo_x(:,1) = agent_velo(:,1);
    all_velo_y(:,1) = agent_velo(:,2);

    % first timestep
    [all_gyrations, traj_init] = calc_all_gyrations(n_agent, n_traj, agent_coor,...
        agent_velo, sigma, box_length, repul_strength, friction, noise,...
        timestep, n_vsteps, v_repul_type);

    [new_coor, new_velo] = next_timestep( agent_coor, ...
        agent_coor, agent_velo, timestep, ...
        n_agent, all_gyrations, traj_init, n_traj, sigma, box_length, repul_strength, ...
        friction, noise, repul_type, synthetic);
    everything_coor_x(:,2) = mod(new_coor(:,1),box_length);
    everything_coor_y(:,2) = mod(new_coor(:,2), box_length);
    all_velo_x(:,2) = new_velo(:,1);
    all_velo_y(:,2) = new_velo(:,2);

    %next timesteps
    for step=3:n_steps
        agent_coor = [everything_coor_x(:,step-1) everything_coor_y(:,step-1)];
        past_agent_coor = [everything_coor_x(:,step-2) everything_coor_y(:,step-2)];
        agent_velo = [all_velo_x(:,step-1) all_velo_y(:,step-1)];
        [all_gyrations, traj_init] = calc_all_gyrations(n_agent, n_traj, agent_coor,...
            agent_velo, sigma, box_length, repul_strength, friction, noise,...
            timestep, n_vsteps, v_repul_type);

        [new_coor, new_velo] = next_timestep( agent_coor, ...
            past_agent_coor, agent_velo, timestep, ...
            n_agent, all_gyrations, traj_init, n_traj, sigma, box_length, repul_strength, ...
            friction, noise, repul_type, synthetic);
        everything_coor_x(:,step) = mod(new_coor(:,1), box_length);
        everything_coor_y(:,step) = mod(new_coor(:,2), box_length);
        all_velo_x(:,step) = new_velo(:,1);
        all_velo_y(:,step) = new_velo(:,2);
        %plot_agents(agent_coor, sigma, box_length, force_init, 1)
        % Plot things -----------------------
%         plot_coor = [everything_coor_x(:,step) everything_coor_y(:,step)];
%         plot_velo = [all_velo_x(:,step) all_velo_y(:,step)];
%         plot_agents(plot_coor, sigma, box_length, plot_velo, 1);
        % -----------------------------------
        prm = "Step " + step + " complete.";
        disp(prm)
        % Find total kinetic energy:
        kinetic_energy = sum( 0.5 * (all_velo_x(:,step).^2 + all_velo_y(:,step).^2));
        disp("Total kinetic energy of step " + step + " is "+ kinetic_energy)
    end
    
    if record
%         vobj = VideoWriter('simulation.avi');
%         vobj.FrameRate = 10;
%         vobj.open;
        fig = figure('visible','on');
        for step = 1:n_steps
            plot_coor = [everything_coor_x(:,step) everything_coor_y(:,step)];
            plot_velo = [all_velo_x(:,step) all_velo_y(:,step)];
            force_repulsive = repulsion(plot_coor, sigma, box_length, ...
                repul_strength, repul_type);
            plot_agents(plot_coor, sigma, box_length, force_repulsive, 1);
%             frame = getframe(fig);
%             writeVideo(vobj,frame);
            prm = "Step " + step + " recorded.";
%             disp(prm)
        end
%         vobj.close;
%         vobj.delete;
    end    
    
    if false
        fig = figure('visible','on');
        for step = 1:n_steps
            plot_coor = [everything_coor_x(:,step) everything_coor_y(:,step)];
            plot_velo = [all_velo_x(:,step) all_velo_y(:,step)];
            plot_agents(plot_coor, sigma, box_length, plot_velo, 1);
        end
    end
   
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% -------------------- Sample from Bivariate Distribution -----------------
% -------------------------------------------------------------------------
function bivar_random = bivariate_normal( samples )
    bivar_random = repmat( [0 0], samples, 1) + ...
        randn(samples, 2) * [1 0; 0 1];
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% -------------------- Virtual Trajectory Creation ------------------------
% -------------------------------------------------------------------------
% i is the agent we are investigating

function [traj_coor, traj_velo, ...
    bound_coor, traj_init_force] = virtual_traj(agent_coor, agent_velo, i, ...
    n_virt_steps, diameter, area, strength, friction, noise, dt, repul_type)
    
    % First element is the real one, so updating number of virtual time
    % steps
    n_virt_steps = n_virt_steps + 1;
    
    % Initialize vectors of path and velocities
    traj_coor = zeros(n_virt_steps,2);
    bound_coor = zeros(n_virt_steps,2);
    traj_velo = zeros(n_virt_steps,2);
    grid_coor = agent_coor;
    
    % Generate random numbers for noise for every step
    virt_steps_noise = bivariate_normal(n_virt_steps);
    % Save the first force to compute gyration later
    traj_init_force = virt_steps_noise(1,:);
    virt_steps_noise = noise * virt_steps_noise;
    
    
    % starting the iteration for the first virtual timestep
    traj_coor(1,:) = agent_coor(i,:);
    traj_velo(1,:) = agent_velo(i,:);
    bound_coor(1,:) = agent_coor(i,:);
    
    f_rep = repulsion_agent(agent_coor, i, diameter, area, strength, repul_type);
    f_langevin = -friction * traj_velo(1,:) + virt_steps_noise(1,:);
    f_tot = f_langevin + f_rep;
    traj_coor(2,:) = traj_coor(1,:) + traj_velo(1,:)* dt + ...
        0.5 * f_tot * (dt^2) ;
    bound_coor(2,:) = mod(traj_coor(2,:), area);
    traj_velo(2,:) = traj_velo(1,:) + f_tot* dt ;
    
    % Update the grid;
    grid_coor(i,:) = mod(traj_coor(2,:), area);
    
    for j=3:n_virt_steps
        % find repulsion force for step
        f_rep = repulsion_agent(grid_coor, i, diameter, area, strength, repul_type);
        f_langevin = -friction * traj_velo(j-1,:) + virt_steps_noise(j,:);
        f_tot = f_rep + f_langevin;
        % update velocity and position of virtual timestep
        traj_coor(j,:) = 2 * traj_coor(j-1,:) - traj_coor(j-2,:) + ...
            f_tot * (dt^2);
        
        %traj_velo(j,:) = (traj_coor(j,:) - traj_coor(j-1,:))/dt;
        traj_velo(j,:) = traj_velo(j-1,:) + f_tot * dt;
        % update the grid:
        bound_coor(j,:) = mod(traj_coor(j,:), area);
        grid_coor(i,:) = mod(traj_coor(j,:), area);
    end
    
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
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



%% ------------------------------------------------------------------------
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



%% ------------------------------------------------------------------------
% ------------- Softcore Repulsion ----------------------------------------
% -------------------------------------------------------------------------
function force_rep = repulsion(agent_coordinates, diameter, area,...
    strength, type)
    if type == "soft"
        force_rep = zeros(length(agent_coordinates), 2);
        for i = 1:size(agent_coordinates, 1)
            for j = 1:size(agent_coordinates, 1)
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
    elseif type == "hard"
        force_rep = zeros(length(agent_coordinates), 2);
        for i = 1:size(agent_coordinates,1)
            for j = 1:size(agent_coordinates,1)
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
                    magnitude = strength;
                    force_rep(i,1) = force_rep(i,1) + magnitude * Dx/ ag_dist;
                    force_rep(i,2) = force_rep(i,2) + magnitude * Dy/ ag_dist;
                end
            end
        end
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% ------------- SOftcore Repulsion for One Agent---------------------------
% -------------------------------------------------------------------------
function force_rep = repulsion_agent(agent_coordinates, i, ...
        diameter, area, strength, type)
    if type == "soft"
        force_rep = zeros(1, 2);
        for j = 1:size(agent_coordinates, 1)
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
%             if i ~= j && ag_dist < 2 * diameter
%                 magnitude = strength*(2*diameter - ag_dist);
%                 force_rep(1) = force_rep(1) + magnitude * Dx/ ag_dist;
%                 force_rep(2) = force_rep(2) + magnitude * Dy/ ag_dist;
%             end
            if i ~= j && ag_dist < 2 * diameter
                magnitude = strength*(2*diameter - ag_dist);
                force_rep(1) = force_rep(1) + magnitude * Dx/ ag_dist;
                force_rep(2) = force_rep(2) + magnitude * Dy/ ag_dist;
            end

        end
    elseif type == "hard"
        force_rep = zeros(1, 2);
        for j = 1:size(agent_coordinates,1)
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
                magnitude = strength;
                force_rep(1) = force_rep(1) + magnitude * Dx/ ag_dist;
                force_rep(2) = force_rep(2) + magnitude * Dy/ ag_dist;
            end
        end
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
               


%% ------------------------------------------------------------------------
% -------------------- Calculate All Gyrations ----------------------------
% -------------------------------------------------------------------------
function [all_gyrations,traj_init] = calc_all_gyrations(n_agent, n_traj, agent_coor,...
    agent_velo, sigma, box_length, repul_strength, friction, noise,...
    timestep, n_vsteps, repul_type)

    all_gyrations = zeros(n_agent, n_traj);
    traj_init = zeros(n_agent,n_traj,2);
    lambdas = zeros(1, n_agent);
    parfor agent_no = 1:n_agent
        %disp('Entering agent')
        for tr=1:n_traj

            [virt_coor, virt_velo, ...
                bound, traj_init_force] = virtual_traj(agent_coor, agent_velo,...
                agent_no, n_vsteps, ...
            sigma, box_length, repul_strength, friction, noise, timestep, repul_type);
            lambdas(agent_no) = lambdas(agent_no)+...
                sqrt((virt_coor(1,1)-virt_coor(n_vsteps,1))^2 + ...
                (virt_coor(1,2)-virt_coor(n_vsteps,2))^2);
            %plot_trajectory(virt_coor,box_length, rand(1,3))
            all_gyrations(agent_no, tr) = calc_gyration(virt_coor);
            traj_init(agent_no, tr,:) = traj_init_force;
        end
        lambdas(agent_no) = lambdas(agent_no)/n_traj;
    end
    disp("Lambda is " + mean(lambdas))
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% -------------------- One Timestep calculation ---------------------------
% -------------------------------------------------------------------------
function [new_coor, new_velo] = next_timestep( agent_coor, ...
    past_agent_coor, agent_velo, timestep, ...
    n_agent, all_gyrations, traj_init, n_traj, diameter, area, ...
    strength, friction, noise, repul_type, synthetic)
    new_coor = zeros(n_agent, 2);
    new_velo = zeros(n_agent, 2);
    for agent=1:n_agent
        if ~ismember(agent, synthetic)
            % Calculate langevin force based on gyration
            f_lang_agent = zeros(1,2);
            agent_mean_gyration = mean(all_gyrations(agent,:));
            for j=1:n_traj
                norm_gyr = log(all_gyrations(agent,j)/agent_mean_gyration);
                f_lang_agent = f_lang_agent + norm_gyr * [traj_init(agent,j,1) traj_init(agent,j,2)];
            end
            f_lang_agent = f_lang_agent/n_traj;
            if isnan(norm_gyr)
                f_lang_agent = -friction * agent_velo(agent,:);
            else
                f_lang_agent = -friction * agent_velo(agent,:) + ...
                    noise * f_lang_agent;
            end
            f_rep = repulsion_agent(agent_coor, agent, diameter, ...
                area, strength, repul_type);
    %         disp("Agent " + agent + ": Langevin is " + f_lang_agent)
            f_tot = f_lang_agent + f_rep;

            % Check the boundaries for position
            past = past_agent_coor(agent,:);
            if (agent_coor(agent, 1)-past_agent_coor(agent,1)) > (0.5 * area)
                past = [(past_agent_coor(agent,1) + area) past_agent_coor(agent,2)];
            elseif (agent_coor(agent, 1)-past_agent_coor(agent,1)) < -(0.5 * area)
                past = [(past_agent_coor(agent,1) - area) past_agent_coor(agent,2)];
            end
            if (agent_coor(agent, 2)-past_agent_coor(agent,2)) > (0.5 * area)
                past = [past_agent_coor(agent,1) (past_agent_coor(agent,2) + area)];
            elseif (agent_coor(agent, 2)-past_agent_coor(agent,2)) < -(0.5 * area)
                past = [past_agent_coor(agent,1) (past_agent_coor(agent,2) - area)];
            end   
            new_coor(agent, :) = 2 * agent_coor(agent, :) - ...
                past + f_tot * (timestep^2);
            % Check the boundaries for velocity
            current = agent_coor(agent, :);
            if (new_coor(agent,1) - current(1)) > (0.5 * area)
                current = [(current(1) + area) current(2)];
            elseif (new_coor(agent,1) - current(1)) < - (0.5 * area)
                current = [(current(1) - area) current(2)];
            end
            if (new_coor(agent,2) - current(2)) > (0.5 * area)
                current = [current(1) (current(2) + area)];
            elseif (new_coor(agent,2) - current(2)) < - (0.5 * area)
                current = [current(1) (current(2) - area) ];
            end
            new_velo(agent, :) = (new_coor(agent,:) - current)/ timestep;
        else
            new_coor(agent, :) = agent_coor(agent, :);
            new_velo(agent, :) = agent_velo(agent, :);
        end
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
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



%% ------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% ----------------------- Runge-Kutta 4 for deterministic model -----------
% -------------------------------------------------------------------------
function [r,v,e] = rk4_abp(r0, v0, e0, d2, gamma, q, c, h, n_steps)
   
    r = zeros(n_steps+1, 2);
    v = zeros(n_steps+1, 2);
    e = zeros(n_steps+1, 1);

    r(1,:) = r0;
    v(1,:) = v0;
    e(1) = e0;

    % Define deterministic ODEs of motion of an active brownian agent
    f = {@(v, e) v;
        @(v, e) d2 * e *v - gamma * v;
        @(v, e) q - c * e - d2 * (v(1)^2 + v(2)^2) * e};

    % % Now solve for number of steps
    for i = 1:n_steps
        k1 =  h * [f{1}(v(i,:),e(i)); f{2}(v(i,:),e(i)); f{3}(v(i,:),e(i)) 0];
        k2 =  h * [f{1}( (v(i,:)+k1(2,:)/2 ),(e(i) + k1(3,1)/2)); ...
            f{2}((v(i,:)+k1(2,:)/2 ),(e(i) + k1(3,1)/2));...
            f{3}((v(i,:)+k1(2,:)/2 ),(e(i) + k1(3,1)/2)) 0];
        k3 =  h * [f{1}( (v(i,:)+k2(2,:)/2 ),(e(i) + k2(3,1)/2)); ...
            f{2}((v(i,:)+k2(2,:)/2 ),(e(i) + k2(3,1)/2));...
            f{3}((v(i,:)+k2(2,:)/2 ),(e(i) + k2(3,1)/2)) 0];
        k4 =  h * [f{1}( (v(i,:)+k3(2,:) ),(e(i) + k3(3,1))); ...
            f{2}((v(i,:)+k3(2,:) ),(e(i) + k3(3,1)));...
            f{3}((v(i,:)+k3(2,:) ),(e(i) + k3(3,1))) 0];
        r(i+1,:) = r(i,:) + (k1(1,:) + 2 .* k2(1,:) + 2.* k3(1,:) + k4(1,:))./6;
        v(i+1,:) = v(i,:) + (k1(2,:) + 2 .* k2(2,:) + 2.* k3(2,:) + k4(2,:))./6;
        e(i+1,:) = e(i,:) + (k1(3,1) + 2 .* k2(3,1) + 2.* k3(3,1) + k4(3,1))./6;
    %     pause(5)
    end


end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% ----------------------- Runge-Kutta 4 general ---------------------------
% -------------------------------------------------------------------------
function [x1_new, x2_new, v1_new, v2_new, e_new] = rk4_step(f, x1, x2, v1, v2, e, h, diameter, ...
    area, strength, repul_type, noise, i, grid_coor)


    f_rep = repulsion_agent(grid_coor, i, diameter, area, strength, repul_type);
    k1 = h * [ f{1}(x1, x2, v1, v2, e);
        f{2}(x1, x2, v1, v2, e);
        f{3}(x1, x2, v1, v2, e, f_rep(1));
        f{4}(x1, x2, v1, v2, e, f_rep(2));
        f{5}(x1, x2, v1, v2, e)];
    k2 = h * [ f{1}(x1 + k1(1)/2, x2+ k1(2)/2, v1+ k1(3)/2, v2+ k1(4)/2, e+ k1(5)/2);
        f{2}(x1 + k1(1)/2, x2+ k1(2)/2, v1+ k1(3)/2, v2+ k1(4)/2, e+ k1(5)/2);
        f{3}(x1 + k1(1)/2, x2+ k1(2)/2, v1+ k1(3)/2, v2+ k1(4)/2, e+ k1(5)/2, f_rep(1));
        f{4}(x1 + k1(1)/2, x2+ k1(2)/2, v1+ k1(3)/2, v2+ k1(4)/2, e+ k1(5)/2, f_rep(2));
        f{5}(x1 + k1(1)/2, x2+ k1(2)/2, v1+ k1(3)/2, v2+ k1(4)/2, e+ k1(5)/2)];
    k3 = h * [ f{1}(x1 + k2(1)/2, x2+ k2(2)/2, v1+ k2(3)/2, v2+ k2(4)/2, e+ k2(5)/2);
        f{2}(x1 + k2(1)/2, x2+ k2(2)/2, v1+ k2(3)/2, v2+ k2(4)/2, e+ k2(5)/2);
        f{3}(x1 + k2(1)/2, x2+ k2(2)/2, v1+ k2(3)/2, v2+ k2(4)/2, e+ k2(5)/2, f_rep(1));
        f{4}(x1 + k2(1)/2, x2+ k2(2)/2, v1+ k2(3)/2, v2+ k2(4)/2, e+ k2(5)/2, f_rep(2));
        f{5}(x1 + k2(1)/2, x2+ k2(2)/2, v1+ k2(3)/2, v2+ k2(4)/2, e+ k2(5)/2)];
    k4 = h * [ f{1}(x1 + k3(1), x2+ k3(2), v1+ k3(3), v2+ k3(4), e+ k3(5));
        f{2}(x1 + k3(1), x2+ k3(2), v1+ k3(3), v2+ k3(4), e+ k3(5));
        f{3}(x1 + k3(1), x2+ k3(2), v1+ k3(3), v2+ k3(4), e+ k3(5), f_rep(1));
        f{4}(x1 + k3(1), x2+ k3(2), v1+ k3(3), v2+ k3(4), e+ k3(5), f_rep(2));
        f{5}(x1 + k3(1), x2+ k3(2), v1+ k3(3), v2+ k3(4), e+ k3(5))];
    x1_new = x1 + (k1(1) + 2*k2(1) + 2*k3(1) + k4(1))/6 ;
    x2_new = x2 + (k1(2) + 2*k2(2) + 2*k3(2) + k4(2))/6 ;
    v1_new = v1 + (k1(3) + 2*k2(3) + 2*k3(3) + k4(3))/6 + noise(1) ;
    v2_new = v2 + (k1(4) + 2*k2(4) + 2*k3(4) + k4(4))/6 + noise(2) ;
    e_new = e + (k1(5) + 2*k2(5) + 2*k3(5) + k4(5))/6 ;

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------