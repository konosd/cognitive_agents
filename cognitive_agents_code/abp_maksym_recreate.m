%% Single Particle Simulations and Visualizations
% parameters, number of agents, trajectories, etc.
n_agent = [50];       %number of agents
n_steps = 2e5;       %number of real steps
n_vsteps = 100;
n_traj = 1;        %number of trajectories
sigma = 1;          %diameter
box_length = 80*sigma;    %area explored

h = 0.005; %timestep = 0.001;     % dt timestep
t = [0:h:(n_steps-1)*h];
virt_t = [0:h:(n_vsteps)*h];
friction = 0.45;     %gamma
temperature = 0.3;  %temperature

D = friction*temperature; %0.01; 
%noise = sqrt(2.0*friction*temperature/timestep);

repul_strength = 1;%20.0;
repul_exp = 10.0;
repul_type = "soft";
v_repul_type = "soft";
pi = 4 * atan(1);

% Add synthetic agents - this is where you define whether an agent will be
% able to move or not. If the number of the agent is in the synthetic
% vector then the agent will not move
synthetic = [];

% parameters for the active brownian agens. Additional ones are: gamma(r)
% additive friction, U(r) potential field, or modify q(r) as an intake
% field. Here, q is only constant. Noise is the same as with passive BP.
q0 = [0 0.5 0.8 1.5 10 5];    % energy intake from the environment
food_radius = 1e6;
food_center = [(box_length*sigma*0.5) (box_length*sigma*0.5)];
% q = @(x1, x2) q0 * ( ((x1-food_center(1))^2 + (x2-food_center(2))^2) < food_radius^2 );
d2 = 3.0;   % conversion rate of internal-to-kinetic energy
c = 1.2;    % dissipation of internal energy



% Potential field
a = 0.0;
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
% % 
% % 
% % agent_coor = initialize_agents(n_agent(k), sigma, box_length);
% agent_coor = [box_length/2, box_length/2];
% agent_velo = [0 0];%zeros(n_agent(k),2);%sqrt(2*D*h)*bivariate_normal(n_agent(k));
% % 
% 
% 
% force_init = repulsion(agent_coor, sigma, box_length, repul_strength, repul_type);
% 
% q = @(x1, x2) q0(iq) * ( ((x1-food_center(1))^2 + (x2-food_center(2))^2) < food_radius^2 );


%% ----------- To visualize one virtual trajectory of one agent -----------

% [my_traj_coor, my_traj_velo, bound, traj_init, traj_e] = virtual_traj(agent_coor, agent_velo, 1, ...
%     n_vsteps(l), sigma, box_length, repul_strength, friction, D, h, v_repul_type, d2, a, c, q);
%     
% 
% fig = figure(1);
% 
% plot_agents(agent_coor, sigma, box_length, force_init, 1)
% plot_trajectory(bound, box_length, rand(1,3))
% circle(food_center, food_radius)
% % title('Simple Brownian Motion')
% xlabel('x')
% ylabel('y')
% xlim([35, 45]);
% ylim([35, 45])
%     
% fig_file = strrep(char(string(['Figures/sp_sim_stoch_food_source_h' , num2str(h) ,  '_g' , ...
%     num2str(friction) , 'char_temp' , num2str(temperature) , '_vsteps' , num2str(n_vsteps(l)) ,'_a' , ...
%     num2str(a) , '_d2' , num2str(d2) , '_c' , num2str(c) , '_q0' , num2str(q0(iq))])), '.', '');
% 
% % set(gcf, 'renderer','Painters')
% % saveas(gcf, fig_file, 'epsc')





%% ---- To visualize everything, solve the full problem, chunk below ------
% % % 


for iq = 1:length(q0)
    for l=1:length(n_vsteps)
        for k=1:length(n_agent)


                    % Filling fraction:
            phi = n_agent(k) * pi * sigma^2 / (4* box_length^2);
            disp("Filling fraction is " + phi)


            %% ------------- Initialization--------------------------------------------
            agent_coor = initialize_agents(n_agent(k), sigma, box_length);
            agent_velo = zeros(n_agent(k),2);

            force_init = repulsion(agent_coor, sigma, box_length, repul_strength, repul_type);

            q = @(x1, x2) q0(iq) * ( ((x1-food_center(1)).^2 + (x2-food_center(2)).^2) < food_radius^2 );

            dir_name = strcat("em_wforce" + n_agent(k) + "_phi"+phi+"_vsteps"+n_vsteps(l)+"_ntraj"+n_traj+"_steps"+n_steps+"_q"+q0(iq));
            mkdir(dir_name)
            
            
            [x, y, u, v, e] = abp_em_cnst_de_solver(n_agent, agent_coor, ...
                    n_steps, sigma, box_length, repul_strength, friction, D, h, repul_type, d2, a, c, q);
            
%             plot_agents(agent_coor, sigma, box_length, force_init, 1.0)
%             hold on
%             for paint_traj = 1:size(agent_coor,1)
%                 plot_trajectory([x(paint_traj,:)' y(paint_traj,:)'] , box_length, rand(1,3))
%             end
                
            incr = 10;
            coordat = zeros(n_agent(k) * n_steps / incr, 5);
            for i=1:n_steps/incr
                coordat(((i-1)*n_agent(k)+1):(i*n_agent(k)) , :) = [x(:,i) y(:,i) ...
                    u(:,i) v(:,i) e(:,i)];
            end

            save(strcat(dir_name, "/coor.dat"), 'coordat', "-ascii")
            

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
        force_rep = zeros(size(agent_coordinates,1), 2);
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
        force_rep = zeros(size(agent_coordinates,1), 2);
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
    elseif type == "exponential"
        force_rep = zeros(size(agent_coordinates,1), 2);
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
                if i ~= j 
                    magnitude = strength / (ag_dist)^2 ;
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
    elseif type == "exponential"
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
            if i ~= j 
                magnitude = strength / (ag_dist^2);
                force_rep(1) = force_rep(1) + magnitude * Dx/ ag_dist;
                force_rep(2) = force_rep(2) + magnitude * Dy/ ag_dist;
            end
        end
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% -------------------- Dissipative Force ----------------------------------
% -------------------------------------------------------------------------
function f_dis = dissipative( agent_coordinates, agent_velocities, ...
    diameter, area)
    
    dis_strength = 1.5;
    repulsion_radius = diameter;
    cutoff_radius = area;

    xij = agent_coordinates(:,1) - agent_coordinates(:,1)' ;
    yij = agent_coordinates(:,2) - agent_coordinates(:,2)' ;
    uij = agent_velocities(:,1) - agent_velocities(:,1)' ;
    vij = agent_velocities(:,2) - agent_velocities(:,2)' ;
    
    rij = sqrt( xij.^2 + yij.^2);
    omegaij = (1 - rij / cutoff_radius).^2;
    omegaij(rij > repulsion_radius) = 0;
    
    magnitude = -dis_strength * omegaij * ( xij.*uij + yij.*vij) ./rij;
    
    f_dis_x = magnitude .* xij ./ rij;
    f_dis_y = magnitude .* yij ./ rij;
    f_dis_x(isnan(f_dis_x)) = 0;
    f_dis_y(isnan(f_dis_y)) = 0;
    
    f_dis_x = sum(f_dis_x,2);
    f_dis_y = sum(f_dis_y,2);
    
    f_dis = [f_dis_x f_dis_x];


end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% -------------------- Solver ---------------------------------------------
% -------------------------------------------------------------------------
% 

function [x, y, u, v, e] = abp_simple_solver(n_agent, agent_coor, ...
    n_steps, diameter, area, strength, friction, D, h, repul_type, d2, a, c, q)
   
    
    % Initialize vectors of path and velocities
    x = zeros(n_agent, n_steps);
    y = zeros(n_agent, n_steps);
    u = zeros(n_agent, n_steps);
    v = zeros(n_agent, n_steps);
    e = zeros(n_agent, n_steps);
    
    x(:,1) = agent_coor(:,1);
    y(:,1) = agent_coor(:,2);
    
    grid_coor = agent_coor;
    
    
    % Generate random numbers for noise for every step
    dw = sqrt(2*D*h) * bivariate_normal(n_steps*n_agent);
    % Save the first force to compute gyration later
    
    

% Runge - Kutta for active particles
    for j = 2:n_steps
        f_rep = repulsion(grid_coor, diameter, area, strength, repul_type);
  
        f_det_x = f_rep(:,1) -friction * u(:,j-1) + ...
            d2 * e(:,j-1) .* u(:,j-1) - 0.5*a*(x(:,j-1)-area/2);
        f_det_y = f_rep(:,2) -friction * v(:,j-1) + ...
            d2 * e(:,j-1) .* v(:,j-1) - 0.5*a*(y(:,j-1)-area/2);

        
        e_det = q(x(:,j-1), y(:,j-1)) - c* e(:,j-1) - d2* e(:,j-1).*(sqrt(v(:,j-1).^2 + u(:,j-1).^2).^2);
        drift = [ u(:,j-1)  v(:,j-1) f_det_x f_det_y e_det];
        
        sk = binornd(1,0.5);
        if sk == 0
            sk = -1;
        end
        volatility = [ zeros(n_agent,1)  zeros(n_agent,1) (sqrt(2*D*h) *bivariate_normal(n_agent)- sk*sqrt(h)) zeros(n_agent,1)];
        k1 = h*drift + volatility ;
        
        %update
        grid_coor = mod([x(:,j-1) y(:,j-1)], area) ;
        
        f_rep = repulsion(grid_coor,  diameter, area, strength, repul_type);
        
        f_det_x = f_rep(:,1) -friction * (u(:,j-1) + k1(:,3)) + ...
            d2 * (e(:,j-1)+k1(:,5)) .*(u(:,j-1) + k1(:,3)) - 0.5*a*( (x(:,j-1)+k1(:,1))-area/2);
        f_det_y = f_rep(:,2) -friction * (v(:,j-1) + k1(:,4)) + ...
            d2 *(e(:,j-1)+k1(:,5)) .*(v(:,j-1) + k1(:,4)) - 0.5*a*( (y(:,j-1)+k1(:,2))-area/2);
        
        
        
        
        e_det = q( (x(:,j-1)+k1(:,1)), (y(:,j-1)+k1(:,2))) - c* (e(:,j-1)+k1(:,5)) -...
            d2* (e(:,j-1)+k1(:,5)).*(sqrt((u(:,j-1) + k1(:,3)).^2 + (v(:,j-1)+k1(:,4)).^2).^2);
        
        drift = [ (u(:,j-1) + k1(:,3)) (v(:,j-1) + k1(:,4))  f_det_x  f_det_y  e_det];
        volatility = [ zeros(n_agent,1)  zeros(n_agent,1) (sqrt(2*D*h) *bivariate_normal(n_agent)+ sk*sqrt(h)) zeros(n_agent,1)];
        k2 = h*drift + volatility;
        
        x(:,j) = mod(x(:,j-1) + 0.5*( k1(:,1) + k2(:,1)), area);
        y(:,j) = mod(y(:,j-1) + 0.5*( k1(:,2) + k2(:,2)), area);
        u(:,j) = u(:,j-1) + 0.5*( k1(:,3) + k2(:,3));
        v(:,j) = v(:,j-1) + 0.5*( k1(:,4) + k2(:,4));
        e(:,j) = e(:,j-1) + 0.5*( k1(:,5) + k2(:,5));
%         disp(strcat("Step " + j + " done."))
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% -------------------- Constant Derivative of Energy depot Solver ---------
% -------------------------------------------------------------------------
% 

function [x, y, u, v, e] = abp_constant_de_solver(n_agent, agent_coor, ...
    n_steps, diameter, area, strength, friction, D, h, repul_type, d2, a, c, q)
   
    
    % Initialize vectors of path and velocities
    x = zeros(n_agent, n_steps);
    y = zeros(n_agent, n_steps);
    u = zeros(n_agent, n_steps);
    v = zeros(n_agent, n_steps);
    e = zeros(n_agent, n_steps);
    
    x(:,1) = agent_coor(:,1);
    y(:,1) = agent_coor(:,2);
    
    grid_coor = agent_coor;
    
    
    % Generate random numbers for noise for every step
    dw = sqrt(2*D*h) * bivariate_normal(n_steps*n_agent);
    % Save the first force to compute gyration later
    
    

% Runge - Kutta for active particles
    for j = 2:n_steps
        f_rep = repulsion(grid_coor, diameter, area, strength, repul_type);
  
        f_det_x = f_rep(:,1) -friction * u(:,j-1) + ...
            d2 * e(:,j-1) .* u(:,j-1) - 0.5*a*(x(:,j-1)-area/2);
        f_det_y = f_rep(:,2) -friction * v(:,j-1) + ...
            d2 * e(:,j-1) .* v(:,j-1) - 0.5*a*(y(:,j-1)-area/2);

        
        e_det = zeros(n_agent,1) ;%q(x(:,j-1), y(:,j-1)) - c* e(:,j-1) - d2* e(:,j-1).*(sqrt(v(:,j-1).^2 + u(:,j-1).^2).^2);
        drift = [ u(:,j-1)  v(:,j-1) f_det_x f_det_y e_det];
        
        sk = binornd(1,0.5);
        if sk == 0
            sk = -1;
        end
        volatility = [ zeros(n_agent,1)  zeros(n_agent,1) (sqrt(2*D*h) *bivariate_normal(n_agent)- sk*sqrt(h)) zeros(n_agent,1)];
        k1 = h*drift + volatility ;
        
        %update
        grid_coor = mod([x(:,j-1) y(:,j-1)], area) ;
        
        f_rep = repulsion(grid_coor,  diameter, area, strength, repul_type);
        
        f_det_x = f_rep(:,1) -friction * (u(:,j-1) + k1(:,3)) + ...
            d2 * (e(:,j-1)+k1(:,5)) .*(u(:,j-1) + k1(:,3)) - 0.5*a*( (x(:,j-1)+k1(:,1))-area/2);
        f_det_y = f_rep(:,2) -friction * (v(:,j-1) + k1(:,4)) + ...
            d2 *(e(:,j-1)+k1(:,5)) .*(v(:,j-1) + k1(:,4)) - 0.5*a*( (y(:,j-1)+k1(:,2))-area/2);
        
        
        
        
        e_det = zeros(n_agent,1) ;%q( (x(:,j-1)+k1(:,1)), (y(:,j-1)+k1(:,2))) - c* (e(:,j-1)+k1(:,5)) -...
            %d2* (e(:,j-1)+k1(:,5)).*(sqrt((u(:,j-1) + k1(:,3)).^2 + (v(:,j-1)+k1(:,4)).^2).^2);
        
        drift = [ (u(:,j-1) + k1(:,3)) (v(:,j-1) + k1(:,4))  f_det_x  f_det_y  e_det];
        volatility = [ zeros(n_agent,1)  zeros(n_agent,1) (sqrt(2*D*h) *bivariate_normal(n_agent)+ sk*sqrt(h)) zeros(n_agent,1)];
        k2 = h*drift + volatility;
        
        x(:,j) = mod(x(:,j-1) + 0.5*( k1(:,1) + k2(:,1)), area);
        y(:,j) = mod(y(:,j-1) + 0.5*( k1(:,2) + k2(:,2)), area);
        u(:,j) = u(:,j-1) + 0.5*( k1(:,3) + k2(:,3));
        v(:,j) = v(:,j-1) + 0.5*( k1(:,4) + k2(:,4));
        e(:,j) = q(x(:,j), y(:,j))./(c+d2*(sqrt(v(:,j).^2 + u(:,j).^2).^2));%    e(:,j-1) + 0.5*( k1(:,5) + k2(:,5));
%         disp(strcat("Step " + j + " done."))

        
    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%% ------------------------------------------------------------------------
% -------------------- E.M. Integrator Cnst de/dt Solver ------------------
% -------------------------------------------------------------------------
% 

function [x, y, u, v, e] = abp_em_cnst_de_solver(n_agent, agent_coor, ...
    n_steps, diameter, area, strength, friction, D, h, repul_type, d2, a, c, q)
   
    
    % Initialize vectors of path and velocities
    x = zeros(n_agent, n_steps);
    y = zeros(n_agent, n_steps);
    u = zeros(n_agent, n_steps);
    v = zeros(n_agent, n_steps);
    e = zeros(n_agent, n_steps);
    
    x(:,1) = agent_coor(:,1);
    y(:,1) = agent_coor(:,2);
    
    grid_coor = agent_coor;
    

    
        % Euler-Murayama method
    for j=2:n_steps
        
        % Generate random numbers for noise for every step
        dw = sqrt(2*D*h) * bivariate_normal(n_agent);
    
        % find repulsion force for step
        f_rep = repulsion(grid_coor, diameter, area, strength, repul_type);
        f_dis = dissipative( grid_coor, [u(:,j-1) v(:,j-1)], diameter, area);
        %

        f_det_x = f_rep(:,1) + f_dis(:,1) -friction * u(:,j-1) + ...
            d2 * e(:,j-1).* u(:,j-1) - 0.5*a*(x(:,j-1)-area/2);
        f_det_y = f_rep(:,2) + f_dis(:,2) -friction * v(:,j-1) + ...
            d2 * e(:,j-1).* v(:,j-1) - 0.5*a*(y(:,j-1)-area/2);
        
        %         f2 = -friction *db + virt_steps_noise(j,:);
        %         disp(f_langevin-f2)
        %         pause(1)
        % update velocity and position of virtual timestep
        
        x(:,j) = mod( (x(:,j-1) + h*u(:,j-1)), area);
        y(:,j) = mod( (y(:,j-1) + h*v(:,j-1)), area);
        u(:,j) = u(:,j-1) + h*f_det_x + dw(:,1);
        v(:,j) = v(:,j-1) + h*f_det_y + dw(:,2);
        e(:,j) = q(x(:,j),y(:,j))./(c+d2*(sqrt(v(:,j).^2 + u(:,j).^2).^2));
        
        % Update grid
        grid_coor = [x(:,j) y(:,j)];
    end

end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------








%% ------------------------------------------------------------------------
% ------------- Function to plot screenshot--------------------------------
% -------------------------------------------------------------------------
function plot_agents(agent_coordinates, diameter, box_length, forces, scl)
    theta = 0:0.01:2*pi;
    for k = 1: size(agent_coordinates,1)
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




