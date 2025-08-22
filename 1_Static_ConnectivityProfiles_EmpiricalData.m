% clean up env for clean run
clearvars; close all; clc

%%%%% INPUTS  %%%%%
   
    %Seed Structure: voxel-wise timeseries (e.g., striatum) 
    
    %execute script in terminal to extract the Seed's voxel-wise timeseries and save them as CSVs:

      % #!/bin/tcsh -xef
      % set subjList = (subj1 subj2 subj3 etc.)
      % foreach subj ($subjList)
      % 3dmaskdump -noijk  -mask rStriatum_Mask.nii.gz {$subj}_Clean_rfMRI_REST1.nii.gz > {$subj}_rStriatum_Clean_TimeSeries_.1D
      % 1dcat -csvout {$subj}_rStriatum_Clean_TimeSeries_.1D > {$subj}_rStriatum_Clean_TimeSeries.csv
      % end
      
    %Target Structures: ROI-wise mean timeseries

    %execute script in terminal to extract the Targets' mean timeseries and save them as CSVs:
    
      % #!/bin/tcsh -xef
      % set subjList = (subj1 subj2 subj3 etc.)
      % set ROIList = (1 2 3 etc.)
      % foreach subj ($subjList)
      % foreach ROI ($ROIList)
      % 3dmaskave -quiet -overwrite -mask {$ROI}.nii {$subj}_Clean_rfMRI_REST1.nii.gz > {$subj}_Clean_rfMRI_REST1_{$ROI}_ts.1D
      % end
      % 1dcat -csvout -overwrite {$subj}_*_Clean_TimeSeries.1D > {$subj}_CorticalROIs_Clean_TimeSeries.csv
      % end
    
%%%%% DOWNLOADSs  %%%%%

    % Download the NIFTI and ANALYZE image toolbox for Matlab:
    %(https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)

addpath '/path/to/Downloads/NIfTI_20140122'


%%RERUN FOR EACH GROUP%% 


%% load time series CSV file paths for each subject

P = '/path/to/seed_CSVs';
S = dir(fullfile(P,'*csv*')); 

P2 = '/path/to/target_CSVs';
S2 = dir(fullfile(P2,'*.csv'));


%% I).

%% Set user parameters

Subjects=  ;                   %User sets parameter 
TRs=  ;                        %User sets parameter 
Seed_Voxels=  ;                %User sets parameter
Total_Target_ROIs=  ;          %User sets parameter
%Restricted_Target_ROIs=  ;  %User sets parameter - optional; restrict to only the strongest X connections

%% load the time series data for each subject

%Group 1%
Final_EdgeTimeseries=zeros(Subjects*TRs*Seed_Voxels,Total_Target_ROIs);
%Final_EdgeTimeseries_Top5=zeros(Subjects*TRs*Seed_Voxels,Restricted_Target_ROIs);

%Group 2%
%Final_EdgeTimeseries=zeros(Subjects*TRs*Seed_Voxels,Total_Target_ROIs);
%Final_EdgeTimeseries_Top5=zeros(Subjects*TRs*Seed_Voxels,Restricted_Target_ROIs);

for i = 1:Subjects
    i=i

%load seed voxel-wise timeseries

    F = fullfile(P,S(i).name);
    S(i).data = readtable(F,'NumHeaderLines', 1);
    S(i).data = table2array(S(i).data); 
    S(i).data = double(S(i).data);
    S(i).data = S(i).data';
    S(i).data = S(i).data(any(S(i).data,2),:);   %remove any motion-censored TRs (i.e. rows of all zeros)
    z1 = zscore(S(i).data);                
    S(i).data = [];

%load target ROI mean timeseries

    F2 = fullfile(P2,S2(i).name);
    S2(i).data = readtable(F2,'NumHeaderLines', 1);
    S2(i).data = table2array(S2(i).data); 
    S2(i).data = double(S2(i).data);
    S2(i).data = S2(i).data(any(S2(i).data,2),:);   %remove any motion-censored TRs (i.e. rows of all zeros)
    z2 = zscore(S2(i).data);           
    S2(i).data = [];

%compute the edge time series of each seed voxel-target ROI pair

for x = 1:Seed_Voxels                                  
    for j = 1:Total_Target_ROIs
       S(i).Final_EdgeTimeseries(1:TRs,j) = z1(:,x).*z2(:,j);
    end
end

%compute the time-averaged connectivity profile of each seed voxel for the subject

for g=1:Seed_Voxels 
 AvgConnProf_acrossTime_Group1(g,1:Total_Target_ROIs,i) = atanh(mean(S(i).Final_EdgeTimeseries(g,:,:),3));
 %AvgConnProf_acrossTime_Group2(g,1:Total_Target_ROIs,i) = atanh(mean(S(i).Final_EdgeTimeseries(g,:,:),3));
end
S(i).Final_EdgeTimeseries = [];

end

%Compute the group-averaged time-averaged connectivity profile of each seed voxel

GroupAvg_ConnProf_Group1 = mean(AvgConnProf_acrossTime_Group1,3);
%GroupAvg_ConnProf_Group2 = mean(AvgConnProf_acrossTime_Group2,3);


%%%% Project the connectivity values into an MNI space NIFTI file %%%%

imgRoot = '/path/to/seedmask';
mask = load_untouch_nii(sprintf('%s/Seed_Mask.nii.gz', imgRoot));

in_brain=find(mask.img(:)>0);
full_brain_space=mask.img.*0;
full_brain_space = double(full_brain_space);

for i=1:Total_Target_ROIs

full_brain_space(in_brain)=GroupAvg_ConnProf_Group1(:,i) 
%full_brain_space(in_brain)=GroupAvg_ConnProf_Group2(:,i) 

imgFile = sprintf('%s/Seed_Mask.nii.gz', imgRoot);
gm = load_untouch_nii(imgFile);
new_hdr=gm;
new_hdr.vol=full_brain_space;
new_hdr.img=full_brain_space;

flNm = num2str(i,'Group1_Seed_Voxelwise_FC_with_ROI_00%d');
%flNm = num2str(i,'Group2_Seed_Voxelwise_FC_with_ROI_00%d');

save_untouch_nii(new_hdr,flNm);

end


