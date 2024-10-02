clc
clear all
close all

fid2=fopen('logfile.txt');
while feof(fid2)==0
    for i = 1:5
        temp=fgetl(fid2);
    end
    maxiteration = sscanf(temp,"%*s %*s %*s %*s %d");
    temp=fgetl(fid2);
    temp=fgetl(fid2);
    dumpfrequency = sscanf(temp,"%*s %*s %*s %d");
    temp=fgetl(fid2);
    temp=fgetl(fid2);
    dt = sscanf(temp,"%*s %*s %f");
    break;
end
fclose(fid2);

mkdir(fullfile(pwd,'images'));

fprintf('rendering images\n');

F1=figure('color','w','visible','off');
set(gcf,'Position',[50 50 1000 1000])
hold on;

ax = gca;
set(ax,'box','on');
axis equal;
set(ax,'Fontsize',26,'XMinorTick','on','YMinorTick','on','ZMinorTick','on');
xlabel('x','Fontsize',30);
ylabel('y','Fontsize',30);
zlabel('z','Fontsize',30);
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
view(35,10);

fprintf('time step = ');
kk = 1;
for n = 0:100*dumpfrequency:maxiteration
    filename = sprintf("dump%08d.off",n);
    fid = fopen(filename);

    if (fid < 0)
        break;
    end

    while feof(fid)==0
        temp=fgetl(fid);
        temp=fgetl(fid);
        [nV,nF,nE]=strread(temp, '%d %d %d');
        
        x = zeros(nV,1);
        y = x;
        z = x;
        tri = zeros(nF,3);
        facealphadata = zeros(nV,1);
        edgealphadata = zeros(nV,1);
        
        for i=1:nV
            temp=fgetl(fid);
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
            temp=fgetl(fid);
            tri(i,:)=sscanf(temp, '%*d %d %d %d');
        end
        tri=tri+1;
    end
    fclose(fid);

    pp1 = patch('Faces',tri,'Vertices',[x,y,z],'FaceColor',"#EDB120",'FaceLighting','gouraud',...
        'FaceVertexAlphaData',facealphadata,'AlphaDataMapping','none','FaceAlpha','interp','EdgeAlpha',0.0);
    pp2 = patch('Faces',tri,'Vertices',[x,y,z],'FaceColor',"#EDB120",'FaceLighting','gouraud',...
        'FaceVertexAlphaData',edgealphadata,'AlphaDataMapping','none','FaceAlpha',0.0,'EdgeAlpha','interp');
    tx = annotation('textbox',[0.05 0.89 0.1 0.1],'String',sprintf('t = %.1f',n*dt),...
        'Fontsize',30,'EdgeColor','w','Color','k');

    imgfname = sprintf('img%05d.jpg',kk);
    fullname = fullfile(pwd,'images',imgfname);
    print(fullname,'-djpeg','-r300');

    delete(pp1);
    delete(pp2);
    delete(tx);

    fprintf('%d ',n);
    kk = kk + 1;
end
fprintf('\n');