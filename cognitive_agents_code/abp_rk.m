%%
% parameters, number of agents, trajectories, etc.
n_agent = [100];       %number of agents
n_vsteps = [800];      %number of virtual steps
n_steps = 120;       %number of real steps
n_traj = 20;        %number of trajectories
sigma = 1;          %diameter
box_length = 80*sigma;    %area explored

h = 0.01; %timestep = 0.001;     % dt timestep
t = [0:h:(n_steps-1)*h];
virt_t = [0:h:(n_vsteps)*h];
friction = 1;     %gamma
temperature = 1;  %temperature

D = friction*temperature; %0.01; 
%noise = sqrt(2.0*friction*temperature/timestep);

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
q0 = [0];    % energy intake from the environment
food_radius = 1e20;
food_center = [80*sigma*0.5 80*sigma*0.5];
% q = @(x1, x2) q0 * ( ((x1-food_center(1))^2 + (x2-food_center(2))^2) < food_radius^2 );
d2 = 1.0;   % conversion rate of internal-to-kinetic energy
c = 1.0;    % dissipation of internal energy



% Potential field
a = 0;
U = @(x1,x2) 0.5 * a * (x1^2 + x2^2);




%-------------------------------------------------------------------------
%% ------------------ Initialization -------------------------------------
% k = 1; %index for number of agents
% iq = 1; % index for magnitutde of energy influx
% l = 1; % index for number of virtual steps
% 
% % Filling fraction:
% phi = n_agent(k) * pi * sigma^2 / (4* box_length^2);
% disp("Filling fraction is " + phi)
% 
% 
% agent_coor = initialize_agents(n_agent(k), sigma, box_length);
% agent_velo = sqrt(2*D*h)*bivariate_normal(n_agent(k));
% 
% % To start from a previous step
% % agent_coor = [all_x(:,1000) all_y(: , 1000)];
% % agent_velo = [vel_x(:,1000) vel_y(: , 1000)];
% 
% force_init = repulsion(agent_coor, sigma, box_length, repul_strength, repul_type);
% 
% q = @(x1, x2) q0(iq) * ( ((x1-food_center(1))^2 + (x2-food_center(2))^2) < food_radius^2 );


%% ----------- To visualize one virtual trajectory of one agent -----------

% [my_traj_coor, my_traj_velo, bound, traj_init] = virtual_traj(agent_coor, agent_velo, 1, ...
%     n_vsteps(l), sigma, box_length, repul_strength, friction, D, h, v_repul_type, d2, a, c, q);
%     
% 
% fig = figure(1);
% plot_agents(agent_coor, sigma, box_length, force_init, 1)
% plot_trajectory(bound, box_length, rand(1,3))

%% ----------- To visualize all virtual trajectories of one agent ---------
% ----------------------- and
% fig = figure(1);
% plot_agents(agent_coor, sigma, box_length, force_init, 1)
% for tr = 1:n_traj
%    [my_virt_traj, my_virt_velo, bound, traj_init] =  virtual_traj(agent_coor, agent_velo, 1, ...
%      n_vsteps, sigma, box_length, repul_strength, friction, D, h, v_repul_type, d2, a, c, q);
%     plot_trajectory(bound, box_length, rand(1,3))
% end

%% ---------- To visualize all virtual trajectories of all agents ---------
% fig = figure(1);
% plot_agents(agent_coor, sigma, box_length, force_init, 1)
% for agent = 1:n_agent
%     for tr = 1:n_traj
%        [my_virt_traj, my_virt_velo, bound, traj_init] =  virtual_traj(agent_coor, ...
%             agent_velo, agent, ...
%          n_vsteps, sigma, box_length, repul_strength, friction, D,...
%               h, v_repul_type, d2, a, c, q );
%         plot_trajectory(bound, box_length, rand(1,3))
%     end
% end

%% ---- To visualize everything, solve the full problem, chunk below ------
% % 
for iq = 1:length(q0)
    for l=1:length(n_vsteps)
        for k=1:length(n_agent)


                    % Filling fraction:
            phi = n_agent(k) * pi * sigma^2 / (4* box_length^2);
            disp("Filling fraction is " + phi)


            %% ------------- Initialization--------------------------------------------
            agent_coor = initialize_agents(n_agent(k), sigma, box_length);
            agent_velo = zeros(n_agent(k),2);

            % To start from a previous step
            % agent_coor = [all_x(:,1000) all_y(: , 1000)];
            % agent_velo = [vel_x(:,1000) vel_y(: , 1000)];

            force_init = repulsion(agent_coor, sigma, box_length, repul_strength, repul_type);

            q = @(x1, x2) q0(iq) * ( ((x1-food_center(1))^2 + (x2-food_center(2))^2) < food_radius^2 );

            dir_name = strcat("abp_agent" + n_agent(k) + "_phi"+phi+"_vsteps"+n_vsteps(l)+"_ntraj"+n_traj+"_steps"+n_steps+"_q"+q0(iq));
            mkdir(dir_name)

            [all_x, all_y, vel_x, vel_y, lambda] = cef_solver( agent_coor,...
                n_agent(k), n_traj, sigma, ...
                box_length, repul_strength, friction, D, h, n_vsteps(l), ...
                n_steps, repul_type, v_repul_type, false, d2, a, c, q);

            incr = 4;
            coordat = zeros(n_agent(k) * n_steps / incr, 4);
            for i=1:n_steps/4
                coordat(((i-1)*n_agent(k)+1):(i*n_agent(k)) , :) = [all_x(:,i) all_y(:,i) ...
                    vel_x(:,i) vel_y(:,i)];
            end

            save(strcat(dir_name, "/coor.dat"), 'coordat', "-ascii")

            save(strcat(dir_name, "/lambdas.dat"), 'lambda', "-ascii")

        end
    end
end





%% ---------------------- End of main -------------------------------------
% -------------------------------------------------------------------------
%% --------------------- Helper Functions ---------------------------------
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
% ------------- Function to initialize agents------------------------------
% -------------------------------------------------------------------------
function init_agents = initialize_agents(n_agents, agent_diameter, box_dim)
    
    init_agents = zeros(n_agents, 2);

    % --------- First set randomly the agents in the box --------------
    % make sure we have enough space
    theor_max = (fix(box_dim*0.5))^2;     %maximum number of agents
    min_distance = 2.1 * agent_diameter;
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
% -------------------- Virtual Trajectory Creation ------------------------
% -------------------------------------------------------------------------
% i is the agent we are investigating

function [traj_coor, traj_velo, ...
    bound_coor, traj_init_force] = virtual_traj(agent_coor, agent_velo, i, ...
    n_virt_steps, diameter, area, strength, friction, D, h, repul_type, d2, a, c, q)
    
    % First element is the real one, so updating number of virtual time
    % steps
    n_virt_steps = n_virt_steps + 1;
    
    % Initialize vectors of path and velocities
    traj_coor = zeros(n_virt_steps,2);
    bound_coor = zeros(n_virt_steps,2);
    traj_velo = zeros(n_virt_steps,2);
    grid_coor = agent_coor;
    
    % Generate random numbers for noise for every step
    dw = sqrt(2*D*h) * bivariate_normal(n_virt_steps);
    % Save the first force to compute gyration later
    traj_init_force = dw(1,:); %/sqrt(2*D);
    
    
    % starting the iteration for the first virtual timestep
    traj_coor(1,:) = agent_coor(i,:);
    traj_velo(1,:) = agent_velo(i,:);
    bound_coor(1,:) = agent_coor(i,:);
    
%     f_rep = repulsion_agent(agent_coor, i, diameter, area, strength, repul_type)
%     f_langevin = -friction * traj_velo(1,:) + dw(1,:) + ...
%         q(traj_coor(1), traj_coor(2)) * d2 * traj_velo(1,:)/(c+d2*norm(traj_velo(1,:))^2)
%     f_tot = f_langevin + f_rep
%     traj_coor(2,:) = traj_coor(1,:) + traj_velo(1,:)* h + ...
%         0.5 * f_tot * (h^2) ;
%     bound_coor(2,:) = mod(traj_coor(2,:), area);
% %     traj_velo(2,:) = traj_velo(1,:) + f_tot* h ;
% %     db = traj_velo(1,:) + f_tot* dt ;
%     traj_velo(2,:) = (traj_coor(2,:) - traj_coor(1,:))/h;
% %     disp(traj_velo(2,:)-db)
%     
%     % Update the grid;
%     grid_coor(i,:) = mod(traj_coor(2,:), area);
    
% Euler-Murayama method
    for j=2:n_virt_steps
        % find repulsion force for step
        f_rep = repulsion_agent(grid_coor, i, diameter, area, strength, repul_type);
        f_det = f_rep -friction * traj_velo(j-1,:) + ...
        q(traj_coor(1), traj_coor(2)) * d2 * traj_velo(j-1,:)/(c+d2*norm(traj_velo(j-1,:))^2);
%         f2 = -friction *db + virt_steps_noise(j,:);
%         disp(f_langevin-f2)
%         pause(1)
        % update velocity and position of virtual timestep
        traj_velo(j,:) = traj_velo(j-1,:) + h*f_det + ...
            dw(j,:);
        
        traj_coor(j,:) = traj_coor(j-1,:) + h*traj_velo(j-1,:);
%         traj_velo(j,:) = traj_velo(j-1,:) + f_tot * h;
%         disp((abs(traj_velo(j,:) - (traj_velo(j-1,:) + f_tot * dt))>1e-9))
        
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
% -------------------- Calculate All Gyrations ----------------------------
% -------------------------------------------------------------------------
function [all_gyrations,traj_init, lambda] = calc_all_gyrations(n_agent, n_traj, agent_coor,...
    agent_velo, sigma, box_length, repul_strength, friction, D,...
    timestep, n_vsteps, repul_type, d2, a, c, q)

    all_gyrations = zeros(n_agent, n_traj);
    traj_init = zeros(n_agent,n_traj,2);
    lambdas = zeros(1, n_agent);
    parfor agent_no = 1:n_agent
        %disp('Entering agent')
        for tr=1:n_traj

            [virt_coor, virt_velo, ...
                bound, traj_init_force] = virtual_traj(agent_coor, agent_velo,...
                agent_no, n_vsteps, sigma, box_length, repul_strength, ...
                friction, D, timestep, repul_type, d2, a, c, q);
            lambdas(agent_no) = lambdas(agent_no)+...
                sqrt((virt_coor(1,1)-virt_coor(n_vsteps,1))^2 + ...
                (virt_coor(1,2)-virt_coor(n_vsteps,2))^2);
            %plot_trajectory(virt_coor,box_length, rand(1,3))
            all_gyrations(agent_no, tr) = calc_gyration(virt_coor);
            traj_init(agent_no, tr,:) = traj_init_force;
        end
        lambdas(agent_no) = lambdas(agent_no)/n_traj;
    end
    lambda = mean(lambdas);
    disp("Lambda is " + lambda)
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% -------------------- One Timestep calculation ---------------------------
% -------------------------------------------------------------------------
function [new_coor, new_velo] = next_timestep( agent_coor, ...
    past_agent_coor, agent_velo, timestep, n_agent, all_gyrations, ...
    traj_init, n_traj, diameter, area, strength, friction, D, repul_type,...
    d2, a, c, q)
    new_coor = zeros(n_agent, 2);
    new_velo = zeros(n_agent, 2);
    for agent=1:n_agent
        % Calculate langevin force based on gyration
        f_lang_agent = zeros(1,2);
        agent_mean_gyration = mean(all_gyrations(agent,:));
        for j=1:n_traj
            norm_gyr = log(all_gyrations(agent,j)/agent_mean_gyration);
            f_lang_agent = f_lang_agent + norm_gyr * [traj_init(agent,j,1) traj_init(agent,j,2)];
        end
        f_lang_agent = f_lang_agent/n_traj;
        if isnan(norm_gyr)
            f_lang_agent = -friction * agent_velo(agent,:) +...
                q(agent_coor(agent,1), agent_coor(agent,2)) * ...
                d2 * agent_velo(agent,:)/(c+d2*norm(agent_velo(agent,:))^2);
        else
            f_lang_agent = -friction * agent_velo(agent,:) + ...
                sqrt(2*D) * f_lang_agent +...
                q(agent_coor(agent,1), agent_coor(agent,2)) * ...
                d2 * agent_velo(agent,:)/(c+d2*norm(agent_velo(agent,:))^2);
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
        
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% --------------------------- All timesteps -------------------------------
% -------------------------------------------------------------------------
function [everything_coor_x, everything_coor_y, all_velo_x, all_velo_y, lambda] =...
    cef_solver( agent_coor,...
    n_agent, n_traj, sigma, ...
    box_length, repul_strength, friction, D, timestep, n_vsteps, ...
    n_steps, repul_type, v_repul_type, record, d2, a, c, q)
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

    lambda = zeros(n_steps,1);
    
    % first timestep
    [all_gyrations, traj_init, lambda(2)] = calc_all_gyrations(n_agent, n_traj, agent_coor,...
        agent_velo, sigma, box_length, repul_strength, friction, D,...
        timestep, n_vsteps, v_repul_type, d2, a, c, q);

    [new_coor, new_velo] = next_timestep( agent_coor, ...
        agent_coor, agent_velo, timestep, ...
        n_agent, all_gyrations, traj_init, n_traj, sigma, box_length, repul_strength, ...
        friction, D, repul_type, d2, a, c, q);
    everything_coor_x(:,2) = mod(new_coor(:,1),box_length);
    everything_coor_y(:,2) = mod(new_coor(:,2), box_length);
    all_velo_x(:,2) = new_velo(:,1);
    all_velo_y(:,2) = new_velo(:,2);

    %next timesteps
    for step=3:n_steps
        agent_coor = [everything_coor_x(:,step-1) everything_coor_y(:,step-1)];
        past_agent_coor = [everything_coor_x(:,step-2) everything_coor_y(:,step-2)];
        agent_velo = [all_velo_x(:,step-1) all_velo_y(:,step-1)];
        [all_gyrations, traj_init, lambda(step)] = calc_all_gyrations(n_agent, n_traj, agent_coor,...
            agent_velo, sigma, box_length, repul_strength, friction, D,...
            timestep, n_vsteps, v_repul_type, d2, a, c, q);

        [new_coor, new_velo] = next_timestep( agent_coor, ...
            past_agent_coor, agent_velo, timestep, ...
            n_agent, all_gyrations, traj_init, n_traj, sigma, box_length, repul_strength, ...
            friction, D, repul_type, d2, a, c, q);
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
% =--------------------------- Plot Circle --------------------------------
% -------------------------------------------------------------------------
% 
function circle(x,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x(1);
    yunit = r * sin(th) + x(2);
    plot(xunit, yunit);
end


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


