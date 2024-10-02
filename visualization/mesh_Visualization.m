clc
clear all
close all

fid1=fopen('dump00010000.off');

while feof(fid1)==0
    temp=fgetl(fid1);
    temp=fgetl(fid1);
    [nV,nF,nE]=strread(temp, '%d %d %d');
    
    x = zeros(nV,1);
    y = x;
    z = x;
    tri = zeros(nF,3);
    
    for i=1:nV
        temp=fgetl(fid1);
        vertex=sscanf(temp, '%g %g %g');
        x(i)=vertex(1);
        y(i)=vertex(2);
        z(i)=vertex(3);    
    end
    
    for i=1:nF
        temp=fgetl(fid1);
        tri(i,:)=sscanf(temp, '%*d %d %d %d');
    end
    tri=tri+1;

end
fclose(fid1);

fid2=fopen('logfile.txt');
while feof(fid2)==0
    for i = 1:17
        temp=fgetl(fid2);
    end
    if (strcmp(sscanf(temp,"%*s %*s %s"),"outside"))
        particle_position = 1;
    elseif (strcmp(sscanf(temp,"%*s %*s %s"),"inside"))
        particle_position = -1;
    end
    temp=fgetl(fid2);
    R0 = sscanf(temp,"%*s %*s %f, %f, %f");
    temp=fgetl(fid2);
    Rp = sscanf(temp,"%*s %*s %f");

    for i = 1:6
        temp=fgetl(fid2);
    end
    chi = sscanf(temp,"%*s %*s %*s %f");
    if (isempty(chi))
        chi = 0.0;
    end
    break;
end
fclose(fid2);

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
axis([-2.0 2.0 0 2.0 -2.0 2.0]);
view(0,0);

h = 2*chi;
Rpp = Rp - 0.01;
[Rx,Ry,Rz] = sphere(200);
X1 = Rpp*Rx; Y1 = Rpp*Ry; Z1 = Rpp*Rz;
X2 = Rpp*Rx; Y2 = Rpp*Ry; Z2 = Rpp*Rz;
if (particle_position > 0)
    X1(Rz>-1+h) = NaN; Y1(Rz>-1+h) = NaN; Z1(Rz>-1+h) = NaN;
    X2(Rz<=-1+h) = NaN; Y2(Rz<=-1+h) = NaN; Z2(Rz<=-1+h) = NaN;
else
    X1(Rz<1-h) = NaN; Y1(Rz<1-h) = NaN; Z1(Rz<1-h) = NaN;
    X2(Rz>=1-h) = NaN; Y2(Rz>=1-h) = NaN; Z2(Rz>=1-h) = NaN;
end

S1=surf(X1+R0(1),Y1+R0(2),Z1+R0(3));
set(S1,'FaceColor',"#A2142F", ...
  'FaceAlpha',1.0,'FaceLighting','gouraud','EdgeColor','none');
S2=surf(X2+R0(1),Y2+R0(2),Z2+R0(3));
set(S2,'FaceColor',"#0072BD", ...
  'FaceAlpha',1.0,'FaceLighting','gouraud','EdgeColor','none');
pp = patch('Faces',tri,'Vertices',[x,y,z],'FaceColor',"#EDB120",'FaceLighting','gouraud','FaceAlpha',1.0,'EdgeColor','k');
%lightangle(45,30);
title(['R_p = ', num2str(Rp), ' \chi = ', num2str(chi)]);