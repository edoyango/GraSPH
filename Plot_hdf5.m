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
    if dim == 2
        axis([0 75 0 40]);
    elseif dim==3
        axis([0 75 0 14.5 -2 42]);
        view(30,30);
    end

    grid on; box on; axis on;
    %     grid off; box off; axis off;

    index=sprintf('%04d',i);

    %% Loading data
    %Physical particle Data
    rx = h5read(['sph_out',index,'.h5'],'/real/x');
    rv = h5read(['sph_out',index,'.h5'],'/real/v');
    rrho = h5read(['sph_out',index,'.h5'],'/real/rho');
    rp = h5read(['sph_out',index,'.h5'],'/real/p');
    rind = h5read(['sph_out',index,'.h5'],'/real/ind');
    rprocid = h5read(['sph_out',index,'.h5'],'/real/procid');
    ritype = h5read(['sph_out',index,'.h5'],'/real/itype');
    
    % Halo particle data
    hx = h5read(['sph_out',index,'.h5'],'/halo/x');
    hv = h5read(['sph_out',index,'.h5'],'/halo/v');
    hrho = h5read(['sph_out',index,'.h5'],'/halo/rho');
    hp = h5read(['sph_out',index,'.h5'],'/halo/p');
    hind = h5read(['sph_out',index,'.h5'],'/halo/ind');
    hprocid = h5read(['sph_out',index,'.h5'],'/halo/procid');
    hitype = h5read(['sph_out',index,'.h5'],'/halo/itype');

    % Virtual particle data
    vx = h5read(['sph_out',index,'.h5'],'/virt/x');
    vv = h5read(['sph_out',index,'.h5'],'/virt/v');
    vrho = h5read(['sph_out',index,'.h5'],'/virt/rho');
    vp = h5read(['sph_out',index,'.h5'],'/virt/p');
    vind = h5read(['sph_out',index,'.h5'],'/virt/ind');
    vprocid = h5read(['sph_out',index,'.h5'],'/virt/procid');
    vitype = h5read(['sph_out',index,'.h5'],'/virt/itype');

    % Ghost particle data
%     gx = h5read(['sph_out',index,'.h5'],'/ghos/x');
%     gv = h5read(['sph_out',index,'.h5'],'/ghos/v');
%     grho = h5read(['sph_out',index,'.h5'],'/ghos/rho');
%     gp = h5read(['sph_out',index,'.h5'],'/ghos/p');
%     gind = h5read(['sph_out',index,'.h5'],'/ghos/ind');
%     gprocid = h5read(['sph_out',index,'.h5'],'/ghos/procid');
%     gitype = h5read(['sph_out',index,'.h5'],'/ghos/itype');

    %% Choosing which quantities to plot and plotting physical particles
    pquant = sqrt(sum(rv(:,:).^2,1));   % speed
%     pquant = rrho;
%     pquant = rprocid;
%     pquant = [1,0,0]                                      % All particles red
    %         ind = find(a00(2,:)==0); a01 = a01(:,ind); pquant = pquant(ind); a00 = a00(:,ind);
    if dim == 2; b = scatter(rx(1,:),rx(2,:),25,pquant,'.');
    elseif dim == 3; b = scatter3(rx(1,:),rx(2,:),rx(3,:),10,pquant,'.');
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
%     vquant = sqrt(sum(vv(:,:).^2,1));   % speed
    % vquant = a22(1,:);
    % ind = find(vquant>0); vx = vx(:,ind); vquant = vquant(ind);
    if dim==2 ; b = scatter(vx(1,:),vx(2,:),50,vquant,'.');
    elseif dim==3; b = scatter3(vx(1,:),vx(2,:),vx(3,:),25,vquant,'.');
    end

    %% plotting ghost particles
    %         gquant = [0.7 0.7 0.7];
%     gquant = sqrt(a31(4,:).^2+a31(5,:).^2+a31(6,:).^2);   % speed
%     %     gquant = -a33(3,:)/9.81/1600;
%     if dim==2 ; b = scatter(a31(1,:),a31(2,:),100,gquant,'.');
%     elseif dim==3; b = scatter3(a31(1,:),a31(2,:),a31(3,:),25,gquant,'.');
%     end

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
fclose('all')
cd(olddir);