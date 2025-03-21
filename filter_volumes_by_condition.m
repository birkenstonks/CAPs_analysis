% return lists of volumes to retain for loops and nonloops
function [loop_idx, nonloop_idx] = filter_volumes_by_condition(filename, scan_legnth)
    
    % read in tsv filename
    stim_times = readtable(filename, "FileType","text", "Delimiter", "\t");
    % disp(stim_times)
    
    loop_idx = []; % loops, for now
    nonloop_idx = []; % nonloops, for now
    delay = 6; % 6 volumes with a TR of 0.8 --> 4.8 seconds (hemodynamic delay)

    % start and end times of each stimulus
    onsets = stim_times.onset;
    ends = stim_times.onset + stim_times.duration;
    
    % finding the volumes that correspond to a given condition
    time = 0;
    for v=1:scan_legnth
        find_stim = find((onsets < time) & (time < ends));
        try
            condition = stim_times.trial_type{find_stim};
            if strcmp(condition, "loop")
                loop_idx(end+1) = v;
            elseif strcmp(condition, "nonloop")
                nonloop_idx(end+1) = v;
            end
        catch
            x = "hello";
        end
        time = time + 0.8;
    end
    

    % shifting the volumes of interest based on hemodynamic delay
    loop_idx = loop_idx + delay;
    nonloop_idx = nonloop_idx + delay;
    
    % removing volumes outside of scan length (after adding delay/shifting the volumes)
    loop_overflow = loop_idx > scan_legnth;
    nonloop_overflow = nonloop_idx > scan_legnth;

    loop_idx(loop_overflow) = [];
    nonloop_idx(nonloop_overflow) = [];

end