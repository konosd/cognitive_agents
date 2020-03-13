% Analysis
% For analysis
%% Single Particle Simulations and Visualizations
% parameters, number of agents, trajectories, etc.
n_agent = [200];       %number of agents
n_steps = 1e5;       %number of real steps
n_vsteps = [120];%[120, 20, 50, 80];
n_traj = 20;        %number of trajectories
sigma = 1;          %diameter
box_length = 80*sigma;    %area explored

h = 0.1; %timestep = 0.001;     % dt timestep
% t = [0:h*incr:(n_steps-1)*h];
virt_t = [0:h:(n_vsteps)*h];
friction = 1.0;     %gamma
temperature = 1.0;  %temperature

D = friction*temperature; %0.01; 
%noise = sqrt(2.0*friction*temperature/timestep);

thermostat = 8;

repul_strength = 10;%20.0;
repul_exp = 10.0;
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
q0 = [1.3];    % energy intake from the environment
food_radius = 1e6;
food_center = [(box_length*sigma*0.5) (box_length*sigma*0.5)];
% q = @(x1, x2) q0 * ( ((x1-food_center(1))^2 + (x2-food_center(2))^2) < food_radius^2 );
d2 = 1.0;   % conversion rate of internal-to-kinetic energy
c = 1.0;    % dissipation of internal energy


% Potential field
a = 0.0;
U = @(x1,x2) 0.5 * a * (x1^2 + x2^2);



%% ---------------------- Import files ------------------------------------
if false
    my_files = dir('synthabp/n200_vsteps120_ntraj20_steps100000_q1.3_thermostat*');
    filenames = strings(length(my_files),1);

    vel_bins = [-5:0.05:5];
    r_bins = [1:80];
    v_dens = zeros([length(q0) length(vel_bins)]);
    u_dens = zeros([length(q0) length(vel_bins)]);
    
    % Correlations
%     c_t = zeros([length(q0), n_steps/incr] );
%     c_r = zeros([length(q0), length(r_bins)]);
    
    
    c_r = zeros(length(filenames), length(r_bins));
    for j=1:length(filenames)
            filenames(j) = my_files(j).name;

            n_agent_j = n_agent( floor((j-1)/length(n_vsteps) + 1));

            coordat=importdata(strcat( 'synthabp/' + filenames(j) + '/coor.dat'));
%             coordat=importdata(strcat( filenames(j) + '/coor.dat'));
            
            incr = n_steps/(size(coordat,1)/n_agent_j);
            x = zeros(n_agent_j, n_steps/incr);
            y = zeros(n_agent_j, n_steps/incr);
            u = zeros(n_agent_j, n_steps/incr);
            v = zeros(n_agent_j, n_steps/incr);
            t = [0:h*incr:(n_steps-1)*h];
%             e = zeros(n_agent_j, n_steps/incr);

            for i=1:n_steps/incr
                x(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 1);
                y(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 2);
                u(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 3);
                v(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 4);
%                 e(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 5);

            end
%             v_dens(j,:) = [mean( reshape(v(:,100:1000),1,[])' < vel_bins)];
%             v_dens(j,:) = [v_dens(j,1) (v_dens(j,2:length(vel_bins)) - v_dens(j,(1:(length(vel_bins)-1))))];
%             u_dens(j,:) = [mean( reshape(u(:,100:1000),1,[])' < vel_bins)];
%             u_dens(j,:) = [u_dens(j,1) (u_dens(j,2:length(vel_bins)) - u_dens(j,(1:(length(vel_bins)-1))))];
            
%             disp([q0(j) min(min(u)) max(max(u)) min(min(v)) max(max(v))])
            
            c_t(j,:) = time_autocor(u, v);
            
            
            tm = 0.9*n_steps/incr:10:n_steps/incr;
            rdf = zeros(length(tm),1);
            for ind = 1:length(tm)
                rdf(ind) =  mean(rad_dist_function([(x(:, tm(ind))) ...
                    (y(:,tm(ind)))], sigma, box_length,200)');
                c_r(j,:) = c_r(j,:) + agent_autocor(x(:,tm(ind)), y(:,tm(ind)), ...
                    u(:,tm(ind)), v(:,tm(ind)),length(r_bins))';
            end
%             rdf(j,:) = rdf(j,:)/(0.1*n_steps/incr );
%             c_r(j,:) = c_r(j,:)/(0.1*n_steps/incr );
            
            phit(j,:) = polar_order(u, v, n_agent_j);
            sum(phit(j,:))


    end
%     clear x y u v
end

if true
plot(t, movmean(c_t, [0 10]))
hold on
xlabel('Time')
ylabel('Mean velocity autocorrelation')
grid on
end


%% ------------------- Plot histograms of velocity distribution------------

if false
% % Make u velocity density distribution
vel_markers = ["-o", "-+", "-*", "-^", "-v", "-x"];
for iq = 1:length(q0)
    plot(vel_bins, sqrt(v_dens(iq,:).^2 + u_dens(iq,:).^2), vel_markers(iq), 'color', rand(1,3))
    hold on
end
legend(strcat("q0 ="+ q0))
xlabel('Velocity')
ylabel('p(V)')
grid on
end

%% ----------------------Analysis for autocorrelations etc.----------------

% 
if false
    
    % Looping through the files
    vel_markers = ["-o", "-+", "-*", "-^"];%, "-v", "-x"];
    color_labels = rand(length(filenames),3);
    ax1 = nexttile;
    for j = 1:length(filenames)
        
        n_a = n_agent( floor((j-1)/length(n_vsteps) + 1));
        
        
        hold(ax1, 'on')
%         set(gca, 'YScale', 'log')
        plot(ax1, t(2:end), movmean(c_t(j,2:end), [100 0]), 'color', color_labels(j,:))
        grid on
    %     subplot(2,1,1)   
    end
    title('Active Time velocity auto-correlation')
    ylabel('Time Velocity autocorrelation')
    xlim([h*4+h n_steps*h])
    xlabel('Time t [s]')
    legend(strcat("n ="+ n_agent))
    hold(ax1,'off')
    
    ax2 = nexttile;
    for j = 1:length(filenames)
        
        hold(ax2,'on')
        set(gca, 'XScale', 'log')
        plot(ax2, r_bins, smoothdata(c_r(j,:)),vel_markers(j), 'color',  color_labels(j,:))
        grid on
    %     subplot(2,1,2, 'XScale','log')   
    end
    title('Active Agent velocity auto-correlation')
    ylabel('Agent Velocity autocorrelation')
    xlim([r_bins(1) r_bins(length(r_bins))])
    xlabel('Distance r [m]')
    legend(strcat("n ="+ n_agent))
    hold(ax2,'off')
        
end

% 

%% -----------------------Radial Distribution Function --------------------

% first all n_agents for one q
if false
    dr = [0:0.1:19.9];
    fig = figure;
    fig.Units = 'centimeters';
    fig.Position(3) = 20;
    fig.Position(4) = 6;
    
    
    for i=1:length(n_agent)
        plot(dr, (movmean(rdf(i,:), [20 0])), 'linewidth', 1.4)
        hold on
    end 
    legend(strcat("n_{a} =" + n_agent))
    set(gca, 'fontsize',12)
    set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02))
    fig.PaperPositionMode = 'auto';
    print('Figures/trial.eps', '-depsc2', '-r600')
    
end

%% -----------------------Polar Order Parameter  --------------------------

% first all n_agents for one q
if false
    
%     fig = figure;
%     fig.Units = 'centimeters';
%     fig.Position(3) = 20;
%     fig.Position(4) = 6;
%     
%     
%     
%     legend(strcat("n_{a} =" + n_agent))
%     set(gca, 'fontsize',12)
%     set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02))
%     fig.PaperPositionMode = 'auto';
%     print('Figures/trial.eps', '-depsc2', '-r600')
%     
    phi=mean(phit,2)';
    [n,ind] = sort(n_agent);
    plot(n, phi(ind))
end

% -------------------------End of Program ---------------------------------


% Script to prepare plot for publication
% fig = figure;
% fig.Units = 'centimeters';
% fig.Position(3) = 20;
% fig.Position(4) = 6;
% % set(gca, 'fontsize',14, 'linewidth', 1.8)
% print -depsc2 Figures/trial
% set(gca, 'fontsize',12)
% set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02))
% fig.PaperPositionMode = 'auto';
% print('Figures/trial.eps', '-depsc2', '-r600')

%% ------------------ Import files for 2 synthetic active agents ----------
if false 

thermostat = [10, 20, 4, 6, 8];
therm_string = ['t2_', 't4_', 't6_', 't8_', 't10_', 't20_'];
q_string = ['q0_', 'q0.3', 'q0.5', 'q0.7', 'q0.9', 'q1.1'];
repeats = 10;
% r10 = zeros(1000,1);
% r20 = zeros(1000,1);
% r4 = zeros(1000,1);
% r6 = zeros(1000,1);
% r8 = zeros(1000,1);
% r2 = zeros(1000,1);

my_files = dir('synthabp/n200_vsteps120_ntraj20_steps100000_q1.3_thermostat*');
filenames = strings(length(my_files),1);

% for j = 1:length(filenames)    
%     filenames(j) = my_files(j).name;
% %    
%     %thermostat( floor((j-1)/repeats +1  )
%     n_agent_j = n_agent(1);
%        coordat=importdata(strcat( 'synthabp/' + filenames(j) + '/coor.dat'));
%             
%         incr = n_steps/(size(coordat,1)/n_agent_j);
%         x = zeros(n_agent_j, n_steps/incr);
%         y = zeros(n_agent_j, n_steps/incr);
%         u = zeros(n_agent_j, n_steps/incr);
%         v = zeros(n_agent_j, n_steps/incr);
%         t = [0:h*incr:(n_steps-1)*h];
% 
%         for i=1:n_steps/incr
%             x(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 1);
%             y(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 2);
%             u(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 3);
%             v(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 4);
% %                 e(:,i) = coordat(((i-1)*n_agent_j+1):(i*n_agent_j) , 5);
% 
%         end
%         r2 = r2 + sqrt( (x(1,:)-x(2,:)).^2 + (y(1,:)-y(2,:)).^2 )'; 
% end
% r2 = r2/10;

% for therm_ix = 1:len(thermo
% 
% 
% end
end
% Now to plot the results of this simulation
if false
    fig = figure;
   markers = ["-o", "-+", "-*", "-^", "-v", "-x"];
%    plot(t,r2, t, r4, t, r6, t, r8, t, r10, t, r20)
   hold on
   plot(t(1:100:end),r2(1:100:end),"-o", t(1:100:end), r4(1:100:end),"-+",...
       t(1:100:end), r6(1:100:end),"-*", t(1:100:end), r8(1:100:end),"-^",...
       t(1:100:end), r10(1:100:end),"-v", t(1:100:end), r20(1:100:end),"-x")
   grid on
   xlabel('Timestep');
   ylabel('Distance');
   hold on
%    legend({'\theta = 2', '\theta = 4', '\theta = 6', '\theta = 8', '\theta = 10', '\theta = 20'});
   legend(strcat("\theta = "+ (num2cell(string([2, 4, 6, 8, 10, 20])))) );
   set(gca, 'fontsize',12)
   fig.Units               = 'centimeters';
    fig.Position(3)         = 20;
    fig.Position(4)         = 8;
    
end

if false
    fig = figure;
    rm = mean([r2 r4 r6 r8 r10 r20]);
    plot( [2 4 6 8 10 20], rm, 'linewidth',2)
    grid on
    xlabel('Thermostat \theta');
    ylabel('\langle d \rangle _{t}');
    set (gca, 'fontsize',12)
    fig.Units = 'centimeters';
    fig.Position(3) = 15;
    fig.Position(4) = 8;
    ylim([0, 35])
end


%% ------------------------- END OF PROGRAM -------------------------------
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
    xij = all_x - all_x';
    
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

function rdf = rad_dist_function(agent_coordinates, sigma, box_length, nbins)
    n = size(agent_coordinates,1);
    extra_bound_coor = [ (agent_coordinates(agent_coordinates(:,1)>(box_length-20),:) - [box_length 0]) ;
            (agent_coordinates(agent_coordinates(:,1)<20,:) + [box_length 0]) ;
            (agent_coordinates(agent_coordinates(:,2)<20,:) + [0 box_length]) ;
            (agent_coordinates(agent_coordinates(:,2)>(box_length-20),:) - [0 box_length]) ];
        enriched_coordinates = [agent_coordinates ; extra_bound_coor];
        
    xij = agent_coordinates(:,1) - enriched_coordinates(:,1)' ;
    yij = agent_coordinates(:,2) - enriched_coordinates(:,2)' ;
    rij = sqrt(xij.^2 + yij.^2);
    
%     rho = size(agent_coordinates,1)*sigma^2*pi()/(4* (box_length)^2);
    rdf = zeros(nbins,1);
    dr = linspace(0,20, nbins+1);
    for i=1:200
        delta = length(rij(rij>dr(i) & rij < (dr(i)+0.1) ))/n;
        delr = pi() * 2 * dr(i) * (0.1);
        rdf(i) = delta* box_length^2 / n /delr;
        rdf(1) = 0;
    end

end


function phit = polar_order(u,v, n_agent)
    nom = u.*mean(u) + v.*mean(v);
    denom = sqrt(u.^2 + v.^2).*sqrt(mean(u).^2 + mean(v).^2);
    phit = sum( nom./denom)/n_agent;
    phit(isnan(phit)) = 0;
end
