%% Visualization

%% Get parameters from running the previous script

np_all = [101];       %number of agents
tau = [100];      %number of virtual steps

nst=n_steps;
type = repul_type;
% Choose to save video and image
save_video=1;
save_figure=0;   

scale=100;  %Scale adjusts arrow size ('30' is an OK one or '0' for no arrows)

nvt=n_traj; %No of virtual trajectories 
box = sigma * 80;
Ar= repul_strength;  %Force strength   
Br = repul_exp;

%% -------------------     Set up visualization     -----------------------

my_files = dir('synthetic_agent_new101_phi0.012395_vsteps100_ntraj360_steps2000_q*');
filenames = strings(length(my_files),1);
for i=1:length(my_files)
   filenames(i) = my_files(i).name;
end


% Video Creation
for j=1:size(filenames)
    % Import data
    dir_in_loop = filenames(j);
    np = np_all( floor((j-1)/length(tau)) + 1);
    coord=importdata(strcat(dir_in_loop + '/coor.dat'));
%     cog_force = importdata(strcat(dir_in_loop + '/cf.dat'));
    % Setup Video Writer
    if save_video
        video_path = char(strcat('/project/home19/kn119/Documents/cognitive_agents/cognitive_agents_code/' +dir_in_loop + '/coor.avi'));
        vobj = VideoWriter(video_path);
        vobj.FrameRate = 10;
        vobj.open;
    end
    fig = figure(1);
    for i = 1:length(coord)/np
        rows = linspace(i*np-np+1,i*np,np);
        vx=coord(rows,3)*scale;
        vy=coord(rows,4)*scale;
%         cfx = cog_force(rows,1)*scale;
%         cfy = cog_force(rows,2)*scale;

        plot(coord(rows,1),coord(rows,2),'.','MarkerSize',30,'color','r');
        hold on
        quiver(coord(rows,1),coord(rows,2),vx,vy,'AutoScale','off','color','k', ...
            'ShowArrowHead','on','LineWidth',0.5);
%         quiver(coord(rows,1),coord(rows,2),cfx,cfy, 'AutoScale','off','color','blue', ...
%             'ShowArrowHead','on','LineWidth',1);
        line('XData', [0 0], 'YData', [0 box], 'LineStyle', '-', ...
            'LineWidth', 1, 'Color','k');
        line('XData', [box box], 'YData', [0 box], 'LineStyle', '-', ...
            'LineWidth', 1, 'Color','k');
        line('XData', [0 box], 'YData', [box box], 'LineStyle', '-', ...
            'LineWidth', 1, 'Color','k');
        line('XData', [0 box], 'YData', [0 0], 'LineStyle', '-', ...
            'LineWidth', 1, 'Color','k');
        axis equal
        grid on;
        axis([0 box 0 box]);
        title(i)
        drawnow
        clear rows
        hold off

         if save_video
             frame = getframe(fig);
             writeVideo(vobj,frame);
         end
    end

    if save_figure
        saveas(strcat(dir_in_loop,"/coor.png"));
    end
    if save_video
        vobj.close;
        vobj.delete;
    end
end


% Image Creation
% for j=1:size(filenames)
%     % Import data
%     dir_in_loop = filenames(j);
%     np = np_all( floor((j-1)/1) + 1);
%     coord=importdata(strcat(dir_in_loop + '/coor.dat'));
%     fig = figure(1);
%     for i = length(coord)/np:length(coord)/np
%         rows = linspace(i*np-np+1,i*np,np);
% 
%         plot
%         drawnow
%         clear rows
%         hold off
% 
%          if save_video
%              frame = getframe(fig);
%              writeVideo(vobj,frame);
%          end
%     end
% 
%     if save_figure
%         saveas(strcat(dir_in_loop,"/coor.png"));
%     end
%     if save_video
%         vobj.close;
%         vobj.delete;
%     end
% end






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
    grid off;
    xlim([0, box_length]);
    ylim([0, box_length])
    hold off;
end