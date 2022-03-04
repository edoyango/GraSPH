%------------------------------------------------------------------------
% This program was developed to plot the output_data from MCGSPH Program.
% The program supports full color representation and also algorithm for 
% making animation data file [Developed by Ha H. BUI, Monash University]
% Modified for pure binary input [Edward Yang, 2020]

%% START PROGRAM
clear; close all; clc

%% program parameters
vid =0; % 1 = make video, 0 = no vide
input_precision = 'double'; % input file precision 'single' for single, 'double' for double

%Define working folder
path = '~/SPH_basic/outputdata';

%parameters
dim = 3;
ini_or_act = 1;
numprocs = 2048;

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

%% Initial configuration
colours = rand(numprocs,3);
if ini_or_act == 0
    % Initial configuration
    fileo0 = fopen([path,'\ini_ind.dat']);             ao0 = fread(fileo0,[3 inf],'*int');
    fileo1 = fopen([path,'\ini_xv.dat']);              ao1 = fread(fileo1,[2*dim inf],input_precision);
    fileo2 = fopen([path,'\ini_state.dat']);           ao2 = fread(fileo2,[2 inf],input_precision);

    pquant = ao0(2,:);

    % plot
    if dim == 2; b = scatter(ao1(1,:),ao1(2,:),50,pquant,'.');
    elseif dim == 3; b = scatter3(ao1(1,:),ao1(2,:),ao1(3,:),50,pquant,'.');
    end
    nshots = length(dir('f_xv*'));

else

    for i=1:80
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
        axis([0 75 0 5 0 40])
        view(30,30)
        grid on; box on; axis on;
    %     grid off; box off; axis off;
    
        index=sprintf('%04d',i); 
    
       %% Loading data
    %     filet1 = fopen(['time',index,'.dat']);              at1 = fread(filet1,[5 inf],input_precision);
        
    %     %Physical particle Data
        file00 = fopen(['f_ind',index,'.dat']);             a00 = fread(file00,[3 inf],'*int');
        file01 = fopen(['f_xv',index,'.dat']);              a01 = fread(file01,[2*dim inf],input_precision);
%         file02 = fopen(['f_state',index,'.dat']);           a02 = fread(file02,[2 inf],input_precision);
        %Halo Particle Data
        file10 = fopen(['h_ind',index,'.dat']);             a10 = fread(file10,[3 inf],'*int');
        file11 = fopen(['h_xv',index,'.dat']);              a11 = fread(file11,[2*dim inf],input_precision);
%         file12 = fopen(['h_state',index,'.dat']);             a12 = fread(file12,[2 inf],input_precision);
        
    %     %Virtual Particle Data
        file20 = fopen(['v_ind',index,'.dat']);             a20 = fread(file20,[3 inf],'*int');
        file21 = fopen(['v_xv',index,'.dat']);              a21 = fread(file21,[2*dim inf],input_precision);
    %     file22 = fopen(['v_state',index,'.dat']);             a22 = fread(file22,[2 inf],input_precision);
    
        %Ghost particle data6
    %     file30 = fopen(['g_ind',index,'.dat']);             a30 = fread(file20,[3 inf],'*int');
    %     file31 = fopen(['g_xv',index,'.dat']);              a31 = fread(file31,[6 inf],input_precision);
    %     file32 = fopen(['g_rho',index,'.dat']);             a32 = fread(file32,[1 inf],input_precision);   
        
        fclose('all');
        
        %% Choosing which quantities to plot and plotting physical particles
        pquant = sqrt(sum(a01(dim+1:2*dim,:).^2,1));   % speed
%         pquant = a01(3,:);
%         pquant = a00(2,:);
%         pquant = [1,0,0]                                      % All particles red
%         ind = find(a00(2,:)==0); a01 = a01(:,ind); pquant = pquant(ind); a00 = a00(:,ind);
        if dim == 2; b = scatter(a01(1,:),a01(2,:),75,pquant,'.');
        elseif dim == 3; b = scatter3(a01(1,:),a01(2,:),a01(3,:),25,pquant,'.');
        end
        
        hold on;
        
        %% plotting halo particles
%         hquant = a10(3,:);
%         hquant = a11(1,:);
    %     hquant = -a13(3,:); % vertical stress
    %     hquant = a10(3,:);
%         hquant = [0 1 1];
    
%         ind = find(a10(2,:)==0); a11 = a11(:,ind); hquant = hquant(ind); a10 = a10(:,ind);
%         b = scatter(a11(1,:),a11(2,:),150,hquant,'.');
    
        %% plotting virtual particles
        vquant = [0.7 0.7 0.7]; % colour all virutla particles grey
%         vquant = a20(1,:);
    %     vquant = sqrt(a21(4,:).^2+a21(5,:).^2+a21(6,:).^2); % speed
    %     vquant = -a23(3,:); % z-stress
    %     vquant = a22(1,:); % density
    %     vquant = a21(6,:); % z-position
    %     vquant = a21(6,:); % z-velocity
    %     ind = find(a21(3,:)>0); a21 = a21(:,ind);
    %     a21(1,:) = a21(1,:) - 5;
    %     for j = 1:length(a21)
    %         xi = a21(1,j)*cos(theta)-a21(3,j)*sin(theta);
    %         yi = a21(1,j)*sin(theta)+a21(3,j)*cos(theta);
    %         a21(1,j) = xi;
    %         a21(3,j) = yi;
    %     end
%         ind = find(a20(2,:)==1); a21 = a21(:,ind);% vquant = vquant(ind);
%         if dim==2 ; b = scatter(a21(1,:),a21(2,:),100,vquant,'.');
%         elseif dim==3; b = scatter3(a21(1,:),a21(2,:),a21(3,:),100,vquant,'.');
%         end
    
        %% plotting ghost particles
%         gquant = [0.7 0.7 0.7];
    %     gquant = sqrt(a31(4,:).^2+a31(5,:).^2+a31(6,:).^2);   % speed
    %     gquant = -a33(3,:)/9.81/1600;
%         b = scatter(a31(1,:),a31(2,:),100,gquant,'.');
    
    %     caxis([0,1])
    
        % text
    %     t = 2.702230588E-3*i;
    %     btext = text(0,.05, .0375,['\delta/B = ',num2str(t*.01/0.02,'%4.3f')],'fontsize',14,'fontWeight','bold');
    %     title('(c) \mu(I)','FontSize',18,'FontWeight','Bold')
        
        hold off;
        drawnow;
        
        %% same figure as frame to video
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
fclose('all')
cd(olddir);