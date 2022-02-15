%------------------------------------------------------------------------
% This program was developed to plot the output_data from sph Program.
% The program supports full color representation and also algorithm for 
% making animation data file

%% START PROGRAM
clear all; close all; clc


%% program parameters
vid =0; % 1 = make video, 0 = no vide

%Define working folder
path = 'C:\Users\edwar\Documents\outputdata';

%parameters
dim = 3;
ini_or_act = 1;

olddir = cd(path);
%% Creating video object
if vid 
    aviobj = VideoWriter('Tdoor070enonlog','Uncompressed AVI');
    aviobj.FrameRate = 15; open(aviobj);
end

%% Setting figure size and renderer
hf     = figure('Renderer','OpenGL','color','w');
hf.Position = [100 100 1200 700];
axis equal; hold on
% axis labels
% xlabel('x-extent','fontsize',12)
% ylabel('y-extent','fontsize',12)
% zlabel('z-extent','fontsize',12)
% set(gcf,'Color',[1,1,1]);
% set(gca,'Units','Centimeters','Position',[1 1 20 15]); % Full axes figure

%% Plotting
nshots = length(dir('f_xv*'));
if ini_or_act == 0
    %Initial configuration=
    ao1 = load([path,'\ini_xv.dat']);
    ao2 = load([path,'\ini_state.dat']);

    pquant = ao1(5,:); % vx

    if dim == 2
        b = scatter(ao1(:,3),ao1(:,4),50,pquant,'.');
    elseif dim == 3
        b = scatter3(ao1(:,3),ao1(:,4),ao1(:,5),50,pquant,'.');
    end

else
    for i=1:10
        %% formatting plot
        clf(hf,'reset'); axis equal; hold on;

        % colourbar
        cbar = colorbar('WestOutside');
    %     cbar.Position = [0.7 0.85 0.2 0.02];
    %     cbar.Ticks = [0 0.05];
        colormap('jet')
    %     ylabel(cbar,'Pressure (kPa)')
        set(cbar,'FontSize',12)
    %     set(gca,'ColorScale','log')
        set(gca,'fontsize',16)
    
        % axis
        axis([-1 76 -1 6 -1 41])
        view(30,30)
        grid on; box on; axis on;
    %     grid off; box off; axis off;

        index=sprintf('%04d',i); 
        
        %% loading data
        %Physical particle Data
        a01 = load(['f_xv',index,'.dat']);
        a02 = load(['f_state',index,'.dat']);
            
        %Virtual Particle Data
        a21 = load(['v_xv',index,'.dat']);
        % a22 = load(['v_state',index,'.dat']);
        
            %Ghost particle data6
        % a31 = load(['g_xv',index,'.dat']);
        % a32 = load(['g_rho',index,'.dat']);
    
        fclose('all');
        
        %% Choosing which quantities to plot and plotting physical particles
        pquant = a02(:,3);
    
        if dim == 2
            b = scatter(a01(:,3),a01(:,4),25,pquant,'.');
        elseif dim == 3
            b = scatter3(a01(:,3),a01(:,4),a01(:,5),25,pquant,'.');
        end
    
        hold on;
    
        %% plotting virtual particles
        vquant = [0.7 0.7 0.7]; % colour all virutla particles grey
    
        if dim == 2
            b = scatter(a21(:,3),a21(:,4),25,vquant,'.');
        elseif dim == 3
            b = scatter3(a21(:,3),a21(:,4),a21(:,5),25,vquant,'.');
        end
    
        %% plotting ghost particles
%         gquant = [0.7 0.7 0.7];
%     
%         if dim == 2
%             b = scatter(a31(:,3),a31(:,4),25,vquant,'.');
%         elseif dim == 3
%             b = scatter3(a31(:,3),a31(:,4),a31(:,5),25,vquant,'.');
%         end
    
        %% misc formatting
        caxis = [0 1];
    
        % text
        %     t = 2.702230588E-3*i;
        %     btext = text(0,.05, .0375,['\delta/B = ',num2str(t*.01/0.02,'%4.3f')],'fontsize',14,'fontWeight','bold');
        %     title('(c) \mu(I)','FontSize',18,'FontWeight','Bold')

        %% make figure visible, write to video file if needed
        hold off;
        drawnow;

        if vid
                F = getframe(hf);
                writeVideo(aviobj,F);
        end

    end
end

%% close video
if vid
    close(aviobj);
    close(hf);
end
cd(olddir);









