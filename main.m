%% Homework Imaging for Neuroscience

% The goal is to investigate the relationship between brain networks described by 
% considering the glucose metabolic pathway provided by PET and by the functional 
% organization as obtained from resting-state fMRI analysis. 

clearvars
close all
clc

%% Preparing the paths

base_path = pwd;                                     % Your path to the folder
utils_path = fullfile(base_path, 'NIFTI_toolbox');   % Insert your Utils Path here
DATA_path = fullfile(base_path, 'DATASET');          % Insert your path to the dataset here

addpath(genpath(utils_path))
addpath(genpath(DATA_path))

%% Brain masking, TAC extraction and Patlak's graphical method

% getting the list of all subjects' data
files = dir(fullfile(DATA_path, '10*'));
% number of subjects
N_subj = length(files);

figure(1)
figure(2)

for sub = 1:N_subj

    % load the data for the subject sub
    disp(['Working on ', files(sub).name])

    % spm
    
    gm = load_untouch_nii(fullfile(files(sub).folder, files(sub).name, 'T1_2_PET.nii', 'c1T1_2_PET.nii'));
    gm = gm.img;
    gm = double(gm)/255;
    
    wm = load_untouch_nii(fullfile(files(sub).folder, files(sub).name, 'T1_2_PET.nii', 'c2T1_2_PET.nii'));
    wm = wm.img;
    wm = double(wm)/255;
    
    csf = load_untouch_nii(fullfile(files(sub).folder, files(sub).name, 'T1_2_PET.nii', 'c3T1_2_PET.nii'));
    csf = csf.img;
    csf = double(csf)/255;

    % load the motion-corrected dynamic PET image
    dynPET = load_nii(fullfile(files(sub).folder, files(sub).name, 'PET_DYN_MoCo.nii.gz'));
    dynData = dynPET.hdr;
    dynPET = dynPET.img;
    dynPET = double(dynPET);

    % creating the summed PET image
    sumPET = sum(dynPET,4);
    
    %%
    % load atlas
    atlas = load_untouch_nii(fullfile(files(sub).folder, files(sub).name, 'Schaefer2018_100Parcels_7Networks_order_2_PET.nii.gz'));
    atlas = atlas.img;
    atlas = double(atlas);
    
    % 3D binary mask of the brain
    brain_mask = (gm + wm + csf) > 0.1; % voxels more likely to be part of the brain

    % mask of voxels to be excluded from the analysis
    exclusion_mask = (wm>0.9) | (csf>0.9) | (~brain_mask);

    % final mask
    epi_mask = (sumPET > 0) & (~exclusion_mask);
    maskIdx = find(epi_mask);

    %%
    ROIidx = nonzeros(unique(atlas));
    n_ROI = length(ROIidx);
    
    % 4D PET linearization
    [nx,ny,nz,nvol] = size(dynPET);
    PET2D = zeros([nvol, size(maskIdx, 1)]);
    
    for tt = 1:nvol
    
        tmp = squeeze(dynPET(:,:,:,tt));
        PET2D(tt,:) = tmp(epi_mask);
        clear tmp
    
    end
    
    % TACs extraction
    for rr = 1 : n_ROI
    
        % extracting the useful voxels from the ROI of interest
        ROImask = atlas(maskIdx)==rr;
        voxelTACs = PET2D(:,ROImask(:));
    
        ROI_TAC = nanmean(voxelTACs,2);
        % to be consisted with ROI_TAC, we eliminate NaNs from all the TACs
        voxelTACs(isnan(voxelTACs)) = 0;
    
        % filtering non-physiological values
        ROI_TAC(ROI_TAC<0)      = 0;
        ROI_TAC(isnan(ROI_TAC)) = 0;
        % adding a zero at the beginning of the vector
        ROI_TAC = [0; ROI_TAC];
    
        ROIall(rr).Idx    = find(ROImask);
        ROIall(rr).Label  = ROIidx(rr); % the original label, which will be different than the index of the vector, since we removed some ROIs
        ROIall(rr).TAC    = ROI_TAC;
        ROIall(rr).TACs   = [zeros(size(voxelTACs,2),1)';voxelTACs]; % we add zero on top of all the TACs, to be consistent
        
        clear voxelTACs
    end
    
    TACS = [ROIall.TAC];               % to switch from struct to matrix
    all_TACS(:,:,sub) = [ROIall.TAC];  % to store all the tacs for each subject


    %% Patlak
    % loading the vector of PET frame duration [s]
    load('timePET.mat')

    % time grid of dynamic PET
    hemiscan  = cumsum(t) - t/2;
    time = [0 hemiscan];
    
    % loading the image-derived input function
    load IDIF.mat
    
    % visualization of the image-derived input funtion
    figure(1)
    title('Image-Derived Input Function')
    subplot(3,3,sub)
    plot(t_IDIF, IDIF)
    xlabel('Time (s)')
    ylabel('Activity (kBq/mL)')
    
    % the second data-point is below the zero, so we remove the first two
    % datapoints from the IDIF
    IDIF = IDIF(3:end);
    t_IDIF = t_IDIF(3:end);
    
    % Compute the integral of Cp, interpolation to the PET time grid
    Int_Cp = cumtrapz(t_IDIF, IDIF);
    
    % Interpolate on PET time grid and compute X
    Cp_Patlak = interp1(t_IDIF, IDIF, min(max(time, t_IDIF(1)), t_IDIF(end)))';  
    Int_Cp_timePET = interp1(t_IDIF, Int_Cp, min(max(time, t_IDIF(1)), t_IDIF(end)))';
    X_Patlak = Int_Cp_timePET./Cp_Patlak; 

    % compute Y
    figure(2)
    subplot(3,3,sub)
    subtitle(['Patlak - subject ',num2str(sub)])
    hold on
    for ii=1:n_ROI
    
        Y_Patlak(ii,:) = TACS(:,ii)./Cp_Patlak;

        plot(X_Patlak, Y_Patlak(ii,:), '*-');
        xlabel('x\_Patlak')
        ylabel('y\_Patlak')

    end
    plot(X_Patlak,mean(Y_Patlak),'*-k',LineWidth=2) % adding the mean to the plot 
    
    % choosing t_star based on linearity: we select the last 17 points
    idx_tstar = 37;
    
    % compute Ki
    for ii=1:length(ROIidx)
    
        Y_Patlak_ii = Y_Patlak(ii, idx_tstar:end)';
    
        G(:,1) = X_Patlak(idx_tstar:end); 
        G(:,2) = 1;
    
        estimates = inv(G'*G)*G'*Y_Patlak_ii;
        Ki(ii) = estimates(1);
    
    end
    
    clear IDIF t_IDIF time

    %% Ki z-scoring
    Ki_zscore(sub, :) = zscore(Ki);

end % closing the for cycle, working on every subject separately
disp('z-scored ki computed for every subject')
disp(' ')

%% Correlation for Ki values

% correlation across the subjects for each Scheafer parcels
[MC_Ki, ~] = corr(Ki_zscore, 'Type', 'Spearman');

% setting the diagonal (high and not meaningful values) to zero
MC_Ki(eye(n_ROI) == 1) = 0;

% visualization
load labels.mat
figure(3)
subplot(131)
imagesc(MC_Ki)
axis square
colormap turbo
colorbar
subtitle('Ki MC')

%% Time-Activity MC

% initializing the MC_TAC matrix
MC_TAC = zeros(n_ROI, n_ROI, N_subj);

for sub = 1:N_subj
    % computing Euclidean distances between ROIs
    dists = pdist(all_TACS(:,:,sub)');
    % reshaping the distances in a symmetric square matrix
    dists_mat = squareform(dists);

    % scaling
    dists_scaled = (dists_mat-min(dists_mat(:)))/(max(dists_mat(:))-min(dists_mat(:)));
    
    % Similarity: 1 - distance, clip max to avoid atanh(1)
    sim = 1 - dists_scaled;
    sim = min(sim, 0.9999);

    % MC_TAC matrix for subject sub
    MC_TAC(:,:,sub) = atanh(sim); 
end

% averaging the MC_TAC matrices to obtain a group-level MC_TAC matrix
MC_TAC_average = mean(MC_TAC, 3);

% setting the diagonal (high and not meaningful values) to zero
MC_TAC_average(eye(n_ROI) == 1) = 0;

% visualization
figure(3)
subplot(132)
imagesc(MC_TAC_average)
axis square
colormap turbo
colorbar
subtitle('TAC MC')


%% fMRI FC

% initializing 
FC_z = zeros(n_ROI, n_ROI, N_subj);

for sub = 1:N_subj

    % load the FC matrix
    disp(['Extracting the FC matrix of subject ', files(sub).name])
    addpath(genpath(fullfile(files(sub).folder, files(sub).name)));
    load('FC_matrix.mat');

    % applying Fischer z-transform
    FC_z(:,:,sub) = atanh(FC);

end
disp('Done')
disp(' ')

% averaging the FC matrices to obtain a group-level FC matrix
FC_average = mean(FC_z, 3);
% setting the diagonal (high and not meaningful values) to zero
FC_average(eye(n_ROI) == 1) = 0;

% visualization
figure(3)
subplot(133)
imagesc(FC_average)
axis square
colormap turbo
colorbar
subtitle('FC')


%% filtering MC_Ki
sorted_MC_Ki = sort(MC_Ki(:), 'descend');

% Find the threshold corresponding to the top 20%
nElements = numel(MC_Ki);
nTop = round(0.2 * nElements); % 20% of total elements
threshold = sorted_MC_Ki(nTop);

% Create a mask of values greater than or equal to the threshold
mask_MC_Ki = MC_Ki >= threshold;

% Apply the mask to retain only top 20%, set others to zero
MC_Ki = MC_Ki .* mask_MC_Ki;

%% filtering MC_TAC_average
sorted_MC_TAC_average = sort(MC_TAC_average(:), 'descend');

% Find the threshold corresponding to the top 20%
nElements = numel(MC_TAC_average);
nTop = round(0.2 * nElements); % 20% of total elements
threshold = sorted_MC_TAC_average(nTop);

% Create a mask of values greater than or equal to the threshold
mask_MC_TAC_average = MC_TAC_average >= threshold;

% Apply the mask to retain only top 20%, set others to zero
MC_TAC_average = MC_TAC_average .* mask_MC_TAC_average;

%% filtering FC_average
sorted_FC_average = sort(FC_average(:), 'descend');

% Find the threshold corresponding to the top 20%
nElements = numel(FC_average);
nTop = round(0.2 * nElements); % 20% of total elements
threshold = sorted_FC_average(nTop);

% Create a mask of values greater than or equal to the threshold
mask_FC_average = FC_average >= threshold;

% Apply the mask to retain only top 20%, set others to zero
FC_average = FC_average .* mask_FC_average;

%% visualization of the filtered connectivity matrices
figure(4)
subplot(131)
imagesc(MC_Ki)
axis square
colormap jet
colorbar
subtitle('Ki MC')

figure(4)
subplot(132)
imagesc(MC_TAC_average)
axis square
colormap jet
colorbar
subtitle('MC TAC')

figure(4)
subplot(133)
imagesc(FC_average)
axis square
colormap jet
colorbar
subtitle('FC')


%% Node Degree and Weighted Clustering Coefficient

% for MC_Ki
% normalized weights
W_Ki = weight_conversion(weight_conversion(MC_Ki, 'autofix'), 'normalize');
% weighted clustering coefficient
WCC_Ki = clustering_coef_wu(W_Ki);
% node degree
DEG_Ki = sum(W_Ki, 2);

% for MC_TAC_average
% normalized weights
W_TAC = weight_conversion(weight_conversion(MC_TAC_average, 'autofix'), 'normalize');
% weighted clustering coefficient
WCC_TAC = clustering_coef_wu(W_TAC);
% node degree
DEG_TAC = sum(W_TAC, 2);

% for FC_average
% normalized weights
W_FC = weight_conversion(weight_conversion(FC_average, 'autofix'), 'normalize');
% weighted clustering coefficient
WCC_FC = clustering_coef_wu(W_FC);
% node degree
DEG_FC = sum(W_FC, 2);


%% correlation between node degree values of the three matrices 

cor_Ki_TAC = corr(DEG_Ki, DEG_TAC, 'Type', 'Spearman');

cor_Ki_FC = corr(DEG_Ki, DEG_FC, 'Type', 'Spearman');

cor_TAC_FC = corr(DEG_TAC, DEG_FC, 'Type', 'Spearman');

figure(5)
subplot(131)
scatter(DEG_Ki, DEG_TAC)
xlabel('DEG\_Ki')
ylabel('DEG\_TAC')
title('Correlation: MC\_Ki - TAC\_average')
lsline
text(mean(DEG_Ki), max(DEG_TAC), ['\rho = ', num2str(cor_Ki_TAC, '%.2f')])
axis square

figure(5)
subplot(132)
scatter(DEG_Ki, DEG_FC)
xlabel('DEG\_Ki')
ylabel('DEG\_FC')
title('Correlation: MC\_Ki - FC\_average')
lsline
text(mean(DEG_Ki), max(DEG_FC), ['\rho = ', num2str(cor_Ki_FC, '%.2f')])
axis square

figure(5)
subplot(133)
scatter(DEG_TAC, DEG_FC)
xlabel('DEG\_TAC')
ylabel('DEG\_FC')
title('Correlation: TAC\_average - FC\_average')
lsline
text(mean(DEG_TAC), max(DEG_FC), ['\rho = ', num2str(cor_TAC_FC, '%.2f')])
axis square

% cor_TAC_FC is the only (potentially) interesting correlation finding.
% cor_TAC_FC >> cor_Ki_FC --> the graph structure of FC is more related to
% time-series MC than to across-subject MC.


%% Identifying hubs for the three networks and investigating their correlation

% for MC_Ki
sorted_DEG_Ki = sort(DEG_Ki, 'descend');
sorted_WCC_Ki = sort(WCC_Ki, 'descend');
% Find the threshold corresponding to the top 30%
threshold_DEG = sorted_DEG_Ki(round(0.3 * numel(DEG_Ki)));
threshold_WCC = sorted_WCC_Ki(round(0.3 * numel(WCC_Ki)));
% Getting the index of the nodes identified as 'hubs'
hub_Ki = find(DEG_Ki >= threshold_DEG & WCC_Ki >= threshold_WCC);


% for MC_TAC_average
sorted_DEG_TAC = sort(DEG_TAC, 'descend');
sorted_WCC_TAC = sort(WCC_TAC, 'descend');
% Find the threshold corresponding to the top 30%
threshold_DEG = sorted_DEG_TAC(round(0.3 * numel(DEG_TAC)));
threshold_WCC = sorted_WCC_TAC(round(0.3 * numel(WCC_TAC)));
% Getting the index of the nodes identified as 'hubs'
hub_TAC = find(DEG_TAC >= threshold_DEG & WCC_TAC >= threshold_WCC);


% for FC_average
sorted_DEG_FC = sort(DEG_FC, 'descend');
sorted_WCC_FC = sort(WCC_FC, 'descend');
% Find the threshold corresponding to the top 30%
threshold_DEG = sorted_DEG_FC(round(0.3 * numel(DEG_FC)));
threshold_WCC = sorted_WCC_FC(round(0.3 * numel(WCC_FC)));
% Getting the index of the nodes identified as 'hubs'
hub_FC = find(DEG_FC >= threshold_DEG & WCC_FC >= threshold_WCC);

disp(['Hub nodes for MC_Ki: ' num2str(hub_Ki')])
disp(' ')
disp(['Hub nodes for MC_TAC_average: ' num2str(hub_TAC')])
disp(' ')
disp(['Hub nodes for FC_average: ' num2str(hub_FC')])
disp(' ')

disp('Hub nodes for MC_Ki: ')
disp(char(L(hub_Ki)))
disp(' ')
disp('Hub nodes for MC_TAC_average: ')
disp(char(L(hub_TAC)))
disp(' ')
disp('Hub nodes for FC_average: ')
disp(char(L(hub_FC)))
disp(' ')

