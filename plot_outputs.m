clc; clear; close all;


%%%%% Creating a blue/white map for 2D plot of water content
%%%%% distribution %%%%%%
map = zeros(200,3)+1;
map(:,2) = 1;
map(1:50,1) = linspace(1,0,50)';
map(51:end,1) = 0;
map(51:150,2) = linspace(.9,0,100)';
map(151:200,2) = 0;
map(151:200,3) = linspace(1,0.5,50)';


%%%%% Importing data from the input file "inputs.txt" %%%%%%%
fid = fopen('inputs.txt');
data = zeros(5,1);
i = 0;
while 1
    i = i + 1;
    tline = fgetl(fid);
    tline = strsplit(tline,' ');
    if (i == 6)
        break; 
    end
    data(i) = str2double(tline{1});
end
fclose(fid);

Lx = data(1);
Ly = data(2);
Nx = data(3);
Ny = data(4);
tPlot = data(5);


%%%%% import simulated water content from "theta.dat" %%%%%%%
theta = importdata('theta.dat');
theta = reshape(theta,numel(theta(:,1)),Ny,Nx);


% 2D plot of water content distribution and save animation in .avi file %%%

v = VideoWriter('theta_2Dprofile.avi','Motion JPEG AVI');
open(v)

for  n = 1:numel(theta(:,1))
    h=figure(1);
    clf

    imagesc(flipud(squeeze(theta(n,:,:))))
    caxis([0 max(max(max(theta)))])
    colormap(map)
    title(['Time = ',num2str(n)* tPlot / 60.,' min']);
    y=colorbar;
    title(y,'\theta','fontsize',25)
    ylabel('Height [m]')
    xlabel('Horizontal distance [m]')
    set(gca,'fontsize',30)
    set(gca,'ytick',[1 Ny])
    set(gca,'yticklabel',[Ly 0])
    set(gca,'xtick',[1 Nx])
    set(gca,'xticklabel',[0 Lx])    
    

    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
        x0=50;
        y0=50;
        height = 700;
        width = 1000;
        set(gcf,'units','points','position',[x0,y0,width,height])
        frame = getframe(gcf);
        writeVideo(v, frame);
end
close(v)
