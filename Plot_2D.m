%------------------------------------------------------------------------
% This program was developed to plot the output_data from MCGSPH Program.
% The program supports full color representation and also algorithm for 
% making animation data file [Developed by Ha H. BUI, Monash University]

% START PROGRAM
clear all; close all; clc

%Define working folder
path = 'C:\Users\edwar\Documents\outputdata';

%Making movie
% aviobj = VideoWriter('dambreak','MPEG-4');
% aviobj.FrameRate = 20; open(aviobj);

%Setting screen size & color
scrsz  = get(0,'ScreenSize');
hf     = figure;
set(gcf,'Color',[1,1,1]);
set(gca,'Position',[0.05 0.0 0.9 1.0]); % Full axes figure

%Program starts

for i=2
      
% ****** Plot real particle ******
    [index,msg]=sprintf('%04d',i); 
    fname  = [path '\f_xv',index,'.dat']; fname1 = [path '\f_state',index,'.dat']; 
    a = load(fname); a1 = load(fname1);
    
% ****** Plot velocity changed color ******
    b = scatter(a(:,3),a(:,4),25,a1(:,4),'.'); % density
%     b = scatter(a(:,2),a(:,3),25,a1(:,4),'.'); % Pressure   
%     b = scatter(a(:,3),a(:,4),25,sqrt(a(:,5).^2+a(:,6).^2),'.'); % Velocity

    axis equal; hold on; 
%     caxis([970 1030])
    bar = colorbar;
    colormap('jet')

% ****** Plot the vector of particle ******
%     b = quiver(a(:,2),a(:,3),a(:,4),a(:,5),1.0);
%     set(b,'linewidth',1);
%     set(b,{'color'},{'k'})

% ****** Plot virtual boundary only *******
    fname2 = [path '\v_xv',index,'.dat']; a2 = load(fname2); 
    b = scatter(a2(:,3),a2(:,4),5,[0.0 0.0 0.0],'.');
    
% ****** Plot time on screen *******
    time_step = 0.5417607;
    time = time_step*i;
    text(10, 37,['Time elapsed = ',num2str(time,'%6.4f'),'s'],'fontsize',12,'fontWeight','bold')
    
% ****** Format axis *******
%     axis([-2.0 77 -2.0 42]);
    
    drawnow;    
%     F = getframe(hf);
%     writeVideo(aviobj,F);
    hold off;
end
% close(aviobj);
% close(hf);









