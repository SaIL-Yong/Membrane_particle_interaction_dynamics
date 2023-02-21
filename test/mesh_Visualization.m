clc
clear all
close all


fid=fopen('dump00100000.off');

while feof(fid)==0
    temp=fgetl(fid);
    temp=fgetl(fid);
    [nV,nF,nE]=strread(temp, '%d %d %d');
    
    FV.vertices=zeros(nV,3);
    FV.faces=zeros(nF,3);
    
    for i=1:nV
        temp=fgetl(fid);
        vertex=sscanf(temp, '%g %g %g');
        x(i)=vertex(1);
        y(i)=vertex(2);
        z(i)=vertex(3);    
    end
     x=x';
     y=y';
     z=z';
    
    for i=1:nF
        temp=fgetl(fid);
        tri(i,:)=sscanf(temp, '%*d %d %d %d %*d %*d');
    end
    tri=tri+1;

end

F1=figure;
set(gcf,'Position',[50 50 1000 600])
hold on;
axis equal
hold on
pp= patch('Faces',tri,'Vertices',[x,y,z],'FaceColor','y','FaceAlpha',0.8,'EdgeColor','k');
axis equal;
axis([-1.5 1.5 -1.5 1.5 -2.0 1.5]);
view(90,0);
title("dt = 0.01, Rp = 0.31, u = 5.0, rho = 0.1, w/ angle condition");