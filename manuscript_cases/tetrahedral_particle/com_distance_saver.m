function com_distance_saver()

    clc; clear; close all;

    % ------------------------------
    % 1) Read parameters from logfile.txt
    % ------------------------------
    fid2 = fopen('logfile.txt','r');
    if fid2 < 0
        error('Cannot open logfile.txt.');
    end

    % This logic depends on your logfile layout:
    % Skipping lines until we find maxiteration, dumpfrequency, dt
    for i = 1:5
        temp = fgetl(fid2);
    end
    maxiteration = sscanf(temp, "%*s %*s %*s %*s %d");  % read an integer
    
    temp = fgetl(fid2);
    temp = fgetl(fid2);
    dumpfrequency = sscanf(temp, "%*s %*s %*s %d");   % read an integer

    temp = fgetl(fid2);
    temp = fgetl(fid2);
    dt = sscanf(temp, "%*s %*s %f");                  % read a float

    % Possibly skip extra lines in logfile
    for i = 1:8
        temp = fgetl(fid2);
    end
    fclose(fid2);

    % ------------------------------
    % 2) Adjust maxiteration if needed
    %    (Your highest file is particle00796000.off => n=796000)
    % ------------------------------
    highest_step = 796000;  % <--- If you know your largest file index is 796000
    if maxiteration > highest_step
        fprintf('Warning: logfile says maxiteration=%d but highest file is %d.\n', ...
            maxiteration, highest_step);
        fprintf('Limiting maxiteration to %d.\n', highest_step);
        maxiteration = highest_step;
    end

    % ------------------------------
    % 3) Prepare output file
    % ------------------------------
    outfile = 'com_distance.txt';
    fid_out = fopen(outfile, 'wt');
    if fid_out < 0
        error('Could not open "%s" for writing.', outfile);
    end
    fprintf(fid_out, '%% time  distance\n');  % optional header

    % ------------------------------
    % 4) Main loop: compute and save COM distance
    % ------------------------------
    for n = 0 : dumpfrequency : maxiteration

        % ---------------------------------------------------
        % 4.1 Read "particle%08d.off" for the particle
        % ---------------------------------------------------
        filename_particle = sprintf('particle%08d.off', n);
        fid_p = fopen(filename_particle, 'r');
        if fid_p < 0
            warning('Could not open %s, skipping.', filename_particle);
            continue;
        end

        % Skip "OFF"
        fgetl(fid_p);  % e.g. "OFF"
        % Read nV nF nE
        line2 = fgetl(fid_p);
        vals  = sscanf(line2, '%d %d %d');
        if numel(vals) < 2
            warning('Failed to parse nV,nF from "%s". Skipping.', line2);
            fclose(fid_p);
            continue;
        end
        nV_part = vals(1);
        nF_part = vals(2);

        % Read particle vertices
        px = zeros(nV_part,1);
        py = zeros(nV_part,1);
        pz = zeros(nV_part,1);
        for i = 1:nV_part
            coordline = fgetl(fid_p);
            coords = sscanf(coordline, '%f %f %f');
            px(i) = coords(1);
            py(i) = coords(2);
            pz(i) = coords(3);
        end

        % Skip face data
        for i = 1:nF_part
            fgetl(fid_p);
        end
        fclose(fid_p);

        com_particle = [mean(px), mean(py), mean(pz)];

        % ---------------------------------------------------
        % 4.2 Read "dump%08d.off" for the vesicle
        % ---------------------------------------------------
        filename_vesicle = sprintf('dump%08d.off', n);
        fid_v = fopen(filename_vesicle, 'r');
        if fid_v < 0
            warning('Could not open %s, skipping.', filename_vesicle);
            continue;
        end

        % Skip "OFF"
        fgetl(fid_v);
        line2 = fgetl(fid_v);
        vals  = sscanf(line2, '%d %d %d');
        if numel(vals) < 2
            warning('Failed to parse nV,nF from "%s". Skipping.', line2);
            fclose(fid_v);
            continue;
        end
        nV_ves = vals(1);
        nF_ves = vals(2);

        vx = zeros(nV_ves,1);
        vy = zeros(nV_ves,1);
        vz = zeros(nV_ves,1);
        for i = 1:nV_ves
            coordline = fgetl(fid_v);
            coords = sscanf(coordline, '%f %f %f');
            vx(i) = coords(1);
            vy(i) = coords(2);
            vz(i) = coords(3);
        end

        % Skip face data
        for i = 1:nF_ves
            fgetl(fid_v);
        end
        fclose(fid_v);

        com_vesicle = [mean(vx), mean(vy), mean(vz)];

        % ---------------------------------------------------
        % 4.3 Compute distance & write
        % ---------------------------------------------------
        dist_com = norm(com_particle - com_vesicle);
        time_now = n * dt;

        fprintf(fid_out, '%.6f  %.6f\n', time_now, dist_com);
    end

    % ------------------------------
    % 5) Finalize
    % ------------------------------
    fclose(fid_out);
    fprintf('COM distance saved in "%s".\n', outfile);

end
