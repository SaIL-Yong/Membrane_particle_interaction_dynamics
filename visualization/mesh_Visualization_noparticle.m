clc
clear all
close all

frame = 3499000;
fid1=fopen(sprintf('dump%08d.off',frame));

while feof(fid1)==0
    temp=fgetl(fid1);
    temp=fgetl(fid1);
    [nV,nF,nE]=strread(temp, '%d %d %d');
    
    x = zeros(nV,1);
    y = x;
    z = x;
    tri = zeros(nF,3);
    facealphadata = zeros(nV,1);
    edgealphadata = zeros(nV,1);
    
    for i=1:nV
        temp=fgetl(fid1);
        vertex=sscanf(temp, '%g %g %g');
        x(i)=vertex(1);
        y(i)=vertex(2);
        z(i)=vertex(3);
        if y(i)>0
            facealphadata(i)=1.0;
            edgealphadata(i)=0.5;
        else
            facealphadata(i)=0.2;
            edgealphadata(i)=0.1;
        end
    end
    
    for i=1:nF
        temp=fgetl(fid1);
        tri(i,:)=sscanf(temp, '%*d %d %d %d');
    end
    tri=tri+1;

end
fclose(fid1);

F1=figure('color','w');
set(gcf,'Position',[50 50 1000 1000])
hold on;

ax = gca;
set(ax,'box','on');
axis equal;
set(ax,'Fontsize',26,'XMinorTick','on','YMinorTick','on','ZMinorTick','on');
xlabel('x','Fontsize',30);
ylabel('y','Fontsize',30);
zlabel('z','Fontsize',30);
axis([-2.0 2.0 -2.0 2.0 -2.0 2.0]);
view(0,0);

pp1 = patch('Faces',tri,'Vertices',[x,y,z],'FaceColor',"#EDB120",'FaceLighting','gouraud',...
        'FaceVertexAlphaData',facealphadata,'AlphaDataMapping','none','FaceAlpha','interp','EdgeAlpha',0.0);
pp2 = patch('Faces',tri,'Vertices',[x,y,z],'FaceColor',"#EDB120",'FaceLighting','gouraud',...
        'FaceVertexAlphaData',edgealphadata,'AlphaDataMapping','none','FaceAlpha',0.0,'EdgeAlpha','interp');
%lightangle(45,30);
