%% Visualization
clear; clc

m=1
for g=1:m
    
    np=10
    tau=10
    nst=200

    save_video=0
    save_figure=1  
    
    scale=30  
    
    nvt=20; %No of virtual trajectories 
    box=80;
    sigma=1;
    Ar=20;  %Force strength   

    mkdir out_fig;

    for j=1:size(dir('*coor.dat'))
        % Import data
        %coord=importdata(strcat(dir_name,'_coor.dat'));
        coord=importdata(strcat('0000',num2str(j),'_coor.dat'));
        % Setup Video Writer
        name=strcat(dir_name,'_',num2str(j))
        if save_video
            vobj = VideoWriter([ pwd strcat('\out_fig\',name,'.avi')],'Motion JPEG avi');
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
            axis off
            axis([0 box 0 box]);
            title(i*50)
            drawnow
            clear rows
            hold off

             if save_video
                 frame = getframe(fig);
                 writeVideo(vobj,frame);
             end
            pause(0.1);
        end

        if save_figure
            saveas(gcf,[pwd strcat('\out_fig\',name,'.png')]);
            %saveas(gcf,name)
        end
        if save_video
            vobj.close;
            vobj.delete;
        end
    end
    %cd ..
end