%% Visualization

%fig_param=importdata('fig_param.txt');
%[m,n] = size(fig_param)
m=1
for g=1:m
 %   np=fig_param(g,1)
 %   tau=fig_param(g,2)
 %   nst=fig_param(g,3)
    np=n_agent;
    tau=n_vsteps;
    nst=n_steps
    %prompt='np='
    %np=input(prompt);  % No of particles
    %prompt='tau='
    %tau=input(prompt);  
    %prompt='n_steps='
    %nst=input(prompt);

    save_video=1;
    save_figure=0;
    %prompt='save_video='
    %save_video=input(prompt);
    %prompt='save_figure='
    %save_figure=input(prompt);    
    
    scale=10;  %Scale adjusts arrow size ('30' is an OK one or '0' for no arrows)
    %prompt='arrow_scale='
    %scale=input(prompt);
    
    nvt=n_traj; %No of virtual trajectories 
    box = sigma * 80;
    Ar= repul_strength;  %Force strength   

    dir_name=strcat('n',num2str(np),'nvt',num2str(nvt),'tau',num2str(tau),'s',num2str(sigma),'L',num2str(box),'a',num2str(Ar),'hsnsteps',num2str(nst));
    %dir_name=strcat('n',num2str(np),'nvt',num2str(nvt),'tau',num2str(tau),'s',num2str(sigma),'L',num2str(box),'a',num2str(Ar),'hs_1');

    %mkdir(dir_name);
%     dir_name = 'n100nvt20tau20s1L80a90s'
%     cd(dir_name)

    for j=1:size(dir('*coor.dat'))
        % Import data
        %coord=importdata(strcat(dir_name,'_coor.dat'));
        %coord=importdata(strcat('0000',num2str(j),'_coor.dat'));
        coord=importdata('coor.dat');
        % Setup Video Writer
        name=strcat(dir_name,'_',num2str(j))
        if save_video
            vobj = VideoWriter('coor.avi','Motion JPEG avi');
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
            title(i*50)
            drawnow
            clear rows
            hold off

             if save_video
                 frame = getframe(fig);
                 writeVideo(vobj,frame);
             end
        end

        if save_figure
            saveas('coor.png');
            %saveas(gcf,name)
        end
        if save_video
            vobj.close;
            vobj.delete;
        end
    end
end