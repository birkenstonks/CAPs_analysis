%% Loading atlases 
% Loading a brain mask file and the atlas files for identifying seed regions
maskfile = '/home/zachkaras/fmri/analysis/atlases/MNI152_T1_2mm_brain_mask.nii.gz';
mask = niftiread(maskfile); % loads the full 91x109x91 mask

mni_brain_file = '/home/zachkaras/fmri/analysis/atlases/MNI152_T1_2mm_brain.nii.gz';
mni_brain = niftiread(mni_brain_file);
brain_idx = find(mask>0); % ident fies only the voxels of the brain, not the empty space

% Loading and create an empty template image for writing NIfTI files:
nii_template = load_untouch_nii(maskfile); % load the brain mask for an example in the right (91x109x91) space
nii_template.img = zeros(size(nii_template.img)); % replace the mask with all 0s
empty_brain = nii_template.img; % create an empty_brain image in the template's shape

atlas = niftiread('/home/zachkaras/fmri/analysis/atlases/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz'); % let's pretend this is a result or seed region file
atlas_2d_brain = atlas(brain_idx); % reshapes to 2d within the brain

for i=1:400 % iterating through Schaefer atlas to get indices for each parcel
    parcels{i} = find(atlas_2d_brain == i);
end

%% Loading and filtering fMRI data
datapath = '/home/zachkaras/fmri/fmri_model_data/clean';
files = dir(datapath); for i=3:numel(files); filenames{i-2}=files(i).name; end % find the file names to analyze in the current directory

all_task_timecourses = []; % 
all_task_labels = [];

% read in each subject's data
for i=1:2 %length(filenames)
    tic
    % Loading and shaping current participant's data
    disp(['subject ',num2str(i),' of ',num2str(numel(filenames)),':', filenames{i}]);
    
    data = niftiread(filenames{i}); % loading the data
    
    data_2d = reshape(data, numel(mask), size(data,4)); % reshaping it to 2D (voxels x timepoints)
    data_2d_brain = data_2d(brain_idx,:); % narrowing it down to only the voxels in the brain mask
    
    % calculating average timecourse for each seed region
    timecourses = [];
    for ii=1:numel(parcels) % creates a seed ROIs x timepoints matrix of the average time series in each seed ROI
        timecourses(ii,:) = nanmean(data_2d_brain(parcels{ii},:));
    end
    
    scan_length = size(timecourses, 2);

    % z-score data for each subject
    timecourses = zscore(timecourses);
    
    % find tsv associated with the scan
    participant = filenames{i}(1:end-7);
    taskfile = sprintf("/home/zachkaras/fmri/fmri_model_data/bids_formatted_clean/sub-%s/sub-%s_task-coding_events.tsv", participant, participant);

    % filter to frames only during the task
    [loop_vols, nonloop_vols] = filter_volumes_by_condition(taskfile, scan_length);
    [sorted_labels, idx_to_keep] = sort_and_filter_scan(scan_length, loop_vols, nonloop_vols);
    % creating labels
    % filtered_timecourses = timecourses(:, sorted_vols);

    timecourses(:, ~idx_to_keep) = [];

    % concatenating into large 2d matrices 
    all_task_timecourses = cat(2, all_task_timecourses, timecourses);
    all_task_labels = cat(2, all_task_labels, sorted_labels);
    % disp(i)
    % disp(size(filtered_timecourses))
    % disp(size(sorted_labels))


    % fprintf("%d: volume size (%d), label size (%d)\n", i, size(filtered_timecourses), size(sorted_labels))

    toc
    % break
end



%% Run clustering


% filter data based on conditions


% Running CAPs using MATLAB code from TbCAPS:
% 1. CAP_ConsensusClustering (X, K_range, Subsample_type, Subsample_fraction, n_folds, DistType)

% 2. ComputeClusteringQuality (Consensus, K_range)


% 3. Kmax - valid number of kmax chosen
% 4. handles.K - number of clusters
% 5. n_rep - number of repetitions
% 6. percentages of positive/negative voxels to keep

% CAP_find_activity (for seed based CAPs)

% 7. Run_Clustering ()


% 8. computing similarity

function [sorted_labels, idx_to_keep] = sort_and_filter_scan(scan_length, loop_vols, nonloop_vols) 
    loop_labels = repmat({'loop'}, size(loop_vols));  % 'A' for array A
    nonloop_labels = repmat({'nonloop'}, size(nonloop_vols));  % 'B' for array B

    filtered_vols = [loop_vols, nonloop_vols];
    combined_labels = [loop_labels, nonloop_labels];
    [sorted_vols, sort_idx] = sort(filtered_vols);

    sorted_labels = combined_labels(sort_idx);
    idx_to_keep = zeros(1, scan_length);
    idx_to_keep(sorted_vols) = 1;
end






