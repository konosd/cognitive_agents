% Analysis
% For analysis
incr = 50;
%% Single Particle Simulations and Visualizations
% parameters, number of agents, trajectories, etc.
n_agent = [1000];       %number of agents
n_steps = 2e5;       %number of real steps
n_vsteps = [100];
n_traj = 1;        %number of trajectories
sigma = 1;          %diameter
box_length = 80*sigma;    %area explored

h = 0.005; %timestep = 0.001;     % dt timestep
t = [0:h*incr:(n_steps-1)*h];
virt_t = [0:h:(n_vsteps)*h];
friction = 1/0.45;     %gamma
temperature = 0.3;  %temperature

D = friction*temperature; %0.01; 
%noise = sqrt(2.0*friction*temperature/timestep);

repul_strength = 10;%20.0;
repul_exp = 10.0;
repul_type = "exponential";
v_repul_type = "exponential";
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



%% -----------------Import files ------------------------------------------
if true
    my_files = dir('em_wforce1000_phi0.12272_vsteps100_ntraj1_steps200000_q*');
    filenames = strings(length(my_files),1);

    vel_bins = [-11:0.05:11];
    r_bins = [1:80];
    v_dens = zeros([length(q0) length(vel_bins)]);
    u_dens = zeros([length(q0) length(vel_bins)]);
    
    % Correlations
    c_t = zeros([length(q0), n_steps/incr] );
    c_r = zeros([length(q0), length(r_bins)]);

    for j=1:length(filenames)
            filenames(j) = my_files(j).name;

            n_agent_j = 1000;%n_agent( floor((j-1)/length(n_vsteps) + 1));

            coordat=importdata(strcat(filenames(j) + '/coor.dat'));
            x = zeros(n_agent_j, n_steps/incr);
            y = zeros(n_agent_j, n_steps/incr);
            u = zeros(n_agent_j, n_steps/incr);
            v = zeros(n_agent_j, n_steps/incr);
            e = zeros(n_agent_j, n_steps/incr);

            for i=1:n_steps/incr
                x(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 1);
                y(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 2);
                u(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 3);
                v(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 4);
                e(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 5);

            end
            v_dens(j,:) = [ mean( reshape(v(:,500:2000),1,[])' < vel_bins)];
            v_dens(j,:) = [v_dens(j,1) (v_dens(j,2:length(vel_bins)) - v_dens(j,(1:(length(vel_bins)-1))))];
            u_dens(j,:) = [ mean( reshape(u(:,500:2000),1,[])' < vel_bins)];
            u_dens(j,:) = [u_dens(j,1) (u_dens(j,2:length(vel_bins)) - u_dens(j,(1:(length(vel_bins)-1))))];
            
            disp([q0(j) min(min(u)) max(max(u)) min(min(v)) max(max(v))])
            
            c_t(j,:) = time_autocor(u, v);
            c_r(j,:) = agent_autocor(x(:, n_steps/incr), y(:,n_steps/incr), u(:,n_steps/incr), v(:,n_steps/incr),length(r_bins))';
    

    end
end
%% ------------------- Plot histograms of velocity distribution------------

if true
% % Make u velocity density distribution
vel_markers = ["-o", "-+", "-*", "-^", "-v", "-x"];
for iq = 1:length(q0)
    plot(vel_bins, v_dens(iq,:), vel_markers(iq), 'color', rand(1,3))
    hold on
end
legend(strcat("q0 ="+ q0))
xlabel('v Velocity')
ylabel('p(u)')
grid on
end



%% ----------------------Analysis for autocorrelations etc.----------------

% 
if false
    
    % Looping through the files
    vel_markers = ["-o", "-+", "-*", "-^", "-v", "-x"];
    color_labels = rand(length(filenames),3);
    ax1 = nexttile;
    for j = 1:length(filenames)
        
        n_a = 1000;%n_agent_all( floor((j-1)/length(n_vsteps_all) + 1));
        
        
        hold(ax1, 'on')
        plot(ax1, t(2:end), movmean(c_t(j,2:end), [100 0]), 'color', color_labels(j,:))
        grid on
    %     subplot(2,1,1)   
    end
    title('Active Time velocity auto-correlation')
    ylabel('Time Velocity autocorrelation')
    xlim([h*4+h n_steps*h])
    xlabel('Time t [s]')
    legend(strcat("q0 ="+ q0))
    hold(ax1,'off')
    
    ax2 = nexttile;
    for j = 1:length(filenames)
        
        hold(ax2,'on')
        plot(ax2, r_bins, c_r(j,:),vel_markers(j), 'color',  color_labels(j,:))
        grid on
    %     subplot(2,1,2, 'XScale','log')   
    end
    title('Active Agent velocity auto-correlation')
    ylabel('Agent Velocity autocorrelation')
    xlim([r_bins(1) r_bins(length(r_bins))])
    xlabel('Distance r [m]')
    legend(strcat("q0 ="+ q0))
    hold(ax2,'off')
        
end

% 


% -------------------------End of Program ---------------------------------

%% ------------------------- Functions ------------------------------------
    

function c = time_autocor(vel_x, vel_y)
    [agents, timesteps] = size(vel_x);
    c = zeros(timesteps,1);
    for t = 2:timesteps
        c(t) = (1/agents)* sum( (vel_x(:,2).*vel_x(:,t)+vel_y(:,2).*vel_y(:,t))./...
            ( sqrt(vel_x(:,2).^2 + vel_y(:,2).^2) .* sqrt(vel_x(:,t).^2 + vel_y(:,t).^2)));
    end
end

function c = agent_autocor(all_x, all_y, vel_x, vel_y, dist_expl)
    agents= length(vel_x);
    distance = [1:dist_expl];
    c = zeros(length(distance),1);
    for r = 1:dist_expl
        c(r) = 0;
        for i = 1:agents
            for j = 1:agents
                if i ~= j
                    rij = sqrt( (all_x(i)-all_x(j))^2 + (all_y(i)-all_y(j))^2);
                    if rij <= r
                        c(r) = c(r) + (vel_x(i)*vel_x(j) + vel_y(i)*vel_y(j))/...
                            ( sqrt(vel_x(i)^2 + vel_y(i)^2) * sqrt(vel_x(j)^2 + vel_y(j)^2));
                    end
                end
            end
        end
        c(r) = c(r)/(agents*(agents-1));
    end
end
