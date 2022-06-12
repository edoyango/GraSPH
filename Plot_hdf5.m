%------------------------------------------------------------------------
% This program was developed to plot the output_data from MCGSPH Program.
% The program supports full color representation and also algorithm for 
% making animation data file [Developed by Ha H. BUI, Monash University]
% Modified for pure binary input [Edward Yang, 2020]

%% START PROGRAM
clear; close all; clc

%% program parameters
vid =0; % 1 = make video, 0 = no vide
input_precision = 'single'; % input file precision 'single' for single, 'double' for double

%Define working folder
path = 'outputdata';

%parameters
dim = 3;
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

for i=1
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

    if dim == 2
        axis([0 75 0 40]);
    elseif dim==3
%         axis([0 75 0 14.5 -2 42]);
        view(30,30);
    end

    grid on; box on; axis on;
    %     grid off; box off; axis off;

    index=sprintf('%04d',i);

    %% Loading data
    file_path = ['sph_out',index,'.h5'];
    %Physical particle Data
    try
        [rx,rv,rrho,rp,rind,rprocid,rtype] = ...
            read_h5_data(file_path,'real');
    catch
        fprintf('Failed to find real particle data\n');
    end
    
    % Halo particle data
    try
        [vx,vv,vrho,vp,vind,vprocid,vtype] = ...
            read_h5_data(file_path,'virt');
    catch
        fprintf('Failed to find virtual particle data\n');
    end

    % Virtual particle data
    try
        [hx,hv,hrho,hp,hind,hprocid,htype] = ...
            read_h5_data(file_path,'halo');
    catch
        fprintf('Failed to find ghost particle data\n');
    end

    % Ghost particle data
    try
        [gx,gv,grho,gp,gind,gprocid,gtype] = ...
            read_h5_data(['sph_out',index,'.h5'],'ghos');
    catch
        fprintf('Failed to find ghost particle data\n');
    end


    %% Choosing which quantities to plot and plotting physical particles
    pquant = sqrt(sum(rv(:,:).^2,1));   % speed
%     pquant = rprocid;
%     pquant = [1,0,0]                                      % All particles red
%     ind = find(rprocid==0); rx = rx(:,ind); pquant = pquant(ind);
    b = plot_helper(rx,pquant,25,'.');

    hold on;

    %% plotting halo particles
%     hquant = [0 1 1];
%     hquant = hv;

%     ind = find(hprocid==2); hx = hx(:,ind); hquant = hquant(ind); 
%     b = plot_helper(hx,hquant,25,'x');

    %% plotting virtual particles
    vquant = [0.7 0.7 0.7]; % colour all virutla particles gre
%     vquant = vv;
% %     vquant = sqrt(sum(vv(:,:).^2,1));   % speed
%     % vquant = a22(1,:);
%     ind = find(vprocid==2); vx = vx(:,ind);% vquant = vquant(ind);
    b = plot_helper(vx,vquant,25,'^');

    %% plotting ghost particles
%     gquant = [0.7 0.7 0.7];
    gquant = gv;
%     gquant = sqrt(a31(4,:).^2+a31(5,:).^2+a31(6,:).^2);   % speed
%     %     gquant = -a33(3,:)/9.81/1600;
%     ind = find(gprocid==2); gx = gx(:,ind); gquant = gquant(ind);
%     b = plot_helper(gx,gquant,25,'o');

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

%% close video
if vid
    close(aviobj);
    close(hf);
end
fclose('all');
cd(olddir);

%% helper function to read sph particle data
function [x,v,rho,p,ind,procid,type] = read_h5_data(path,group)
    
    x = h5read(path,['/',group,'/x']);
    v = h5read(path,['/',group,'/v']);
    rho = h5read(path,['/',group,'/rho']);
    p = h5read(path,['/',group,'/p']);
    ind = h5read(path,['/',group,'/ind']);
    procid = h5read(path,['/',group,'/procid']);
    type = h5read(path,['/',group,'/type']);
end

function b = plot_helper(x,data,sz,symbol)
    if size(x,1)==2
        b = scatter(x(1,:),x(2,:),sz,data,symbol);
    elseif size(x,1)==3
        b = scatter3(x(1,:),x(2,:),x(3,:),sz,data,symbol);
    end
end