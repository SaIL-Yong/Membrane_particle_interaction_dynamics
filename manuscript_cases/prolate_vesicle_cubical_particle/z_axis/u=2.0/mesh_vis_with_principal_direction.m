% Clear workspace and close figures
clc
clear all
close all

% --------------------------
% 1) Read parameters from logfile
% --------------------------
fid2 = fopen('logfile.txt');
while feof(fid2) == 0
    for i = 1:5
        temp = fgetl(fid2);
    end
    maxiteration = sscanf(temp, "%*s %*s %*s %*s %d");
    temp = fgetl(fid2);
    temp = fgetl(fid2);
    dumpfrequency = sscanf(temp, "%*s %*s %*s %d");
    temp = fgetl(fid2);
    temp = fgetl(fid2);
    dt = sscanf(temp, "%*s %*s %f");
    for i = 1:8
        temp = fgetl(fid2);
    end
    break;
end
fclose(fid2);

% --------------------------
% 2) Read principal axes data
% --------------------------
fid_axes = fopen('transformed_principal_axes.txt', 'r');
fgetl(fid_axes);  % Skip one header line
principal_axes_cell = textscan(fid_axes, '%f %f %f %f %f %f %f %f %f %f');
fclose(fid_axes);
principal_axes_data = cell2mat(principal_axes_cell);

% --------------------------
% 3) Create directory for images
% --------------------------
mkdir(fullfile(pwd, 'images'));
fprintf('rendering images\n');

% --------------------------
% 4) Create figure
% --------------------------
F1 = figure('Color','w','Visible','off');
set(F1, 'Units', 'pixels', 'Position', [50 50 1000 1000], ...
         'PaperPositionMode','auto', ...
         'InvertHardcopy','off');

hold on;
ax = gca;
set(ax, 'box', 'off', 'FontSize', 26, ...
        'XMinorTick','off','YMinorTick','off','ZMinorTick','off');
xlabel('x','FontSize',30); ylabel('y','FontSize',30); zlabel('z','FontSize',30);

axis equal;
axis([-1.40 1.40 -1.40 1.40 -2.80 3.50]);  % your fixed axis limits
ax.XColor = 'none';
ax.YColor = 'none';
ax.ZColor = 'none';

view(120, 10);
lightangle(60, 10);
lightangle(-60, -10);
lighting gouraud;

% Fill figure space & lock axis
set(ax, 'Units','normalized','Position',[0 0 1 1]);
set(ax, 'LooseInset',[0 0 0 0]);
axis manual;

fprintf('time step = ');
kk = 1;

% Overwrite dt if you want
dt = 0.001;

quiverHandles = [];

% ------------------------------------------------------
% Save a "base" camera. We'll do it once, after the figure is created.
% Then, inside the loop, we'll reset to this baseline each time.
% ------------------------------------------------------
drawnow;  % Ensure figure is rendered before we snapshot the camera
baseCamPos   = get(ax,'CameraPosition');
baseCamTgt   = get(ax,'CameraTarget');
baseCamUpVec = get(ax,'CameraUpVector');
baseCamVA    = get(ax,'CameraViewAngle');  % for perspective mode

% If you want perspective:
camproj('perspective');  % We'll also ensure this each time

% --------------------------
% 5) Main loop
% --------------------------
for n = 0:dumpfrequency:maxiteration
    
    % 5.1: Find the current timestep in principal_axes_data
    timestep_idx = find(principal_axes_data(:,1) == n, 1);
    if isempty(timestep_idx)
        continue;
    end

    % 5.2: Read "particle%08d.off"
    filename = sprintf("particle%08d.off", n);
    fid1 = fopen(filename);
    if (fid1 < 0), break; end

    while feof(fid1) == 0
        temp = fgetl(fid1);  % OFF
        temp = fgetl(fid1);  % nV nF nE
        [nV, nF, ~] = strread(temp,'%d %d %d');
        
        x1 = zeros(nV,1);
        y1 = x1; z1 = x1;
        tri1 = zeros(nF,3);

        for i = 1:nV
            temp = fgetl(fid1);
            vertex = sscanf(temp, '%g %g %g');
            x1(i) = vertex(1);  y1(i) = vertex(2);  z1(i) = vertex(3);
        end

        for i = 1:nF
            temp = fgetl(fid1);
            tri1(i,:) = sscanf(temp, '%*d %d %d %d');
        end
        tri1 = tri1 + 1;
    end
    fclose(fid1);

    % 5.3: Delete old principal directions if any
    if ~isempty(quiverHandles)
        delete(quiverHandles);
    end

    % 5.4: Retrieve principal directions
    principal_dir = reshape(principal_axes_data(timestep_idx, 2:end), [3, 3]);
    center_of_mass = [mean(x1), mean(y1), mean(z1)];

    % quiverHandles = [
    %     quiver3(center_of_mass(1), center_of_mass(2), center_of_mass(3), ...
    %             principal_dir(1,1), principal_dir(2,1), principal_dir(3,1), ...
    %             'LineWidth',4, 'MaxHeadSize',0.5, 'Color','r','AutoScale','on','AutoScaleFactor',0.5);
    %     quiver3(center_of_mass(1), center_of_mass(2), center_of_mass(3), ...
    %             principal_dir(1,2), principal_dir(2,2), principal_dir(3,2), ...
    %             'LineWidth',4, 'MaxHeadSize',0.5, 'Color','g','AutoScale','on','AutoScaleFactor',0.5);
    %     quiver3(center_of_mass(1), center_of_mass(2), center_of_mass(3), ...
    %             principal_dir(1,3), principal_dir(2,3), principal_dir(3,3), ...
    %             'LineWidth',4, 'MaxHeadSize',0.5, 'Color','b','AutoScale','on','AutoScaleFactor',0.5);
    % ];

    % 5.5: Read "dump%08d.off"
    filename = sprintf("dump%08d.off", n);
    fid = fopen(filename);
    if (fid < 0), break; end

    while feof(fid) == 0
        temp = fgetl(fid);  % OFF
        temp = fgetl(fid);  % nV nF nE
        [nV, nF, ~] = strread(temp, '%d %d %d');
        
        x = zeros(nV,1);
        y = x; z = x;
        tri = zeros(nF,3);

        for i = 1:nV
            temp = fgetl(fid);
            vertex = sscanf(temp, '%g %g %g');
            x(i) = vertex(1); y(i) = vertex(2); z(i) = vertex(3);
        end

        for i = 1:nF
            temp = fgetl(fid);
            tri(i,:) = sscanf(temp, '%*d %d %d %d');
        end
        tri = tri + 1;
    end
    fclose(fid);

    pp = patch('Faces', tri, 'Vertices', [x,y,z], ...
               'FaceColor',"#C0C0C0", 'FaceLighting','gouraud', ...
               'FaceAlpha',0.2, 'EdgeColor','k','EdgeAlpha',0.3);

    pp1 = patch('Faces', tri1, 'Vertices', [x1,y1,z1], ...
                'FaceColor',"#FFA500", 'FaceLighting','gouraud', ...
                'FaceAlpha',0.6, 'EdgeColor','k','EdgeAlpha',0.5);

    % ------------------------------------------------------
    % 5.6: RESET the camera to the base viewpoint each time
    % ------------------------------------------------------
    set(ax, 'CameraPosition',   baseCamPos);
    set(ax, 'CameraTarget',     baseCamTgt);
    set(ax, 'CameraUpVector',   baseCamUpVec);
    set(ax, 'CameraViewAngle',  baseCamVA);
    camproj('perspective');  % ensure perspective

    % Then apply a consistent zoom factor (e.g. 1.5) from that base
    camzoom(1.5);

    % 5.7: Annotation
    % tx = annotation('textbox',[0.05 0.80 0.1 0.1], ...
    %     'String', sprintf('t=%.1f', n*dt), ...
    %     'FontSize',64, 'EdgeColor','w', 'Color','k');

    % 5.8: Print image
    imgfname = sprintf('img%05d.jpg', kk);
    fullname = fullfile(pwd, 'images', imgfname);
    print(fullname, '-djpeg', '-r400');

    delete(pp);
    delete(pp1);
    % delete(tx);

    fprintf('%d ', n);
    kk = kk + 1;
end

fprintf('\n');
