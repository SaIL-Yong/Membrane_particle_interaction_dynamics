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
    for i = 1:8
        temp=fgetl(fid2);
    end

    break;
end
fclose(fid2);

R_p=0.3;
chi = 0.0;
h = 2*chi;


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
axis([-2.0 2.0 -2.0 2.0 -2.0 2.0]);
view(0,0);

%title(['R_p = ', Rp, ' \chi = ', chi);

fprintf('time step = ');
kk = 1;

filename = sprintf("initialparticle.off");
    fid1 = fopen(filename);

    % if (fid1 < 0)
    %     break;
    % end

    while feof(fid1)==0
        temp=fgetl(fid1);
        temp=fgetl(fid1);
        [nV,nF,nE]=strread(temp, '%d %d %d');
        
        x1 = zeros(nV,1);
        y1 = x1;
        z1 = x1;
        tri1 = zeros(nF,3);

        for i=1:nV
            temp=fgetl(fid1);
            vertex=sscanf(temp, '%g %g %g');
            x1(i)=vertex(1);
            y1(i)=vertex(2);
            z1(i)=vertex(3);    
        end
        
        for i=1:nF
            temp=fgetl(fid1);
            tri1(i,:)=sscanf(temp, '%*d %d %d %d');
        end
        tri1=tri1+1;
    end
fclose(fid1);

for n = 0:dumpfrequency:maxiteration
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

        for i=1:nV
            temp=fgetl(fid);
            vertex=sscanf(temp, '%g %g %g');
            x(i)=vertex(1);
            y(i)=vertex(2);
            z(i)=vertex(3);    
        end
        
        for i=1:nF
            temp=fgetl(fid);
            tri(i,:)=sscanf(temp, '%*d %d %d %d');
        end
        tri=tri+1;
    end
    fclose(fid);

    pp = patch('Faces',tri,'Vertices',[x,y,z],'FaceColor',"	#feff36",'FaceLighting','gouraud','FaceAlpha',0.2,'EdgeColor','k');
    pp1 = patch('Faces',tri1,'Vertices',[x1,y1,z1],'FaceColor',"#f75a22",'FaceLighting','gouraud','FaceAlpha',0.4,'EdgeColor','k');
    tx = annotation('textbox',[0.05 0.89 0.1 0.1],'String',sprintf('t = %.1f',n*dt),...
        'Fontsize',30,'EdgeColor','w','Color','k');

    imgfname = sprintf('img%05d.jpg',kk);
    fullname = fullfile(pwd,'images',imgfname);
    print(fullname,'-djpeg','-r300');

    delete(pp);
    delete(pp1);
    delete(tx);

    fprintf('%d ',n);
    kk = kk + 1;
end
fprintf('\n');