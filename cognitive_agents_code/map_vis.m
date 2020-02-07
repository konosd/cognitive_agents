%% Visualization

%% Get parameters from running the previous script
np=n_agent;
tau=n_vsteps;
nst=n_steps;
type = repul_type;
% Choose to save video and image
save_video=1;
save_figure=0;   

scale=1;  %Scale adjusts arrow size ('30' is an OK one or '0' for no arrows)

nvt=n_traj; %No of virtual trajectories 
box = sigma * 80;
Ar= repul_strength;  %Force strength   
Br = repul_exp;

%% -------------------     Set up visualization     -----------------------



file_name=strcat('new_abp_',num2str(np),'_nvt',num2str(nvt),'_tau',num2str(tau),...
    '_s',num2str(sigma),'_L',num2str(box),'_A',num2str(Ar),'_nst',num2str(nst),...
    '_f', type, '_q',num2str(q0),'_potfield',num2str(a));


for j=1:size(dir('*coor.dat'))
    % Import data
    coord=importdata('coor.dat');
    % Setup Video Writer
    if save_video
        vobj = VideoWriter(strcat(file_name, 'coor.avi'),'MPEG-4');
        vobj.FrameRate = 10;
        vobj.open;
    end
    fig = figure(1);
    for i = 1:length(coord)/np
        rows = linspace(i*np-np+1,i*np,np);
        vx=coord(rows,3)*scale;
        vy=coord(rows,4)*scale;

        plot(coord(rows,1),coord(rows,2),'.','MarkerSize',30,'color','r');
        hold on
        quiver(coord(rows,1),coord(rows,2),vx,vy,'AutoScale','off','color','k', ...
            'ShowArrowHead','on','LineWidth',1);
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
        saveas(strcat(file_name,'coor.png'));
    end
    if save_video
        vobj.close;
        vobj.delete;
    end
end
