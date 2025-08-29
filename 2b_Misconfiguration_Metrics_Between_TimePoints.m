%Uses input from Script 1%

% Script for computing voxel-wise connectivity profile metrics (i.e., aggregate
% divergence, rank order misarrangement, entropy shift) between different timepoints of the same group


clc;
clear;
close all;

softwareRoot = '/Applications/MATLAB_R2020a.app';;
%script requires download of NIFTI and ANALYZE image toolbox for Matlab
%(https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
addpath '/NIfTI_20140122'
addpath(genpath('/MATLAB'))


GroupAvg_Aggregate_Divergence = zeros(Seed_Voxels);          %For storing the final, voxel-wise Aggregate Divergence values
GroupAvg_RankOrderRearrangement = zeros(Seed_Voxels);        %For storing the final, voxel-wise Rank Order Rearrangement values
GroupAvg_Entropy_Difference = zeros(Seed_Voxels);            %For storing the final, voxel-wise Entropy Shift values

  
%%%%%%%%%
%%Compute the change between time points in each of the 3 connectivity profile properties at each voxel in each individual 

for i=1:Subjects

  for b=1:Seed_Voxels
   
%%%Entropy Difference

    h1= abs(AvgConnProf_acrossTime_Group1(b,:,i));
    p1 = bsxfun(@rdivide,h1,sum(h1',1));
    e1 = -nansum(p1.*log2(p1),1);
    Entropy_T1 = e1/log2(Total_Target_ROIs);
      
    
    h2= abs(AvgConnProf_acrossTime_Group2(b,:,i));
    p2 = bsxfun(@rdivide,h2,sum(h2',1));
    e2 = -nansum(p2.*log2(p2),1);
    Entropy_T2 = e2/log2(Total_Target_ROIs);
    
    Subject_Entropy_Difference(b,i) = abs(Entropy_T1-Entropy_T2);
  
%%%Aggregate Divergence (Manhattan Distance)

  for j=1:Total_Target_ROIs
  
    MD(j)=abs(AvgConnProf_acrossTime_Group1(b,j,i)-AvgConnProf_acrossTime_Group2(b,j,i));
  
  end 
  
    MD_total=sum(MD);
    Subject_Aggregate_Divergence(b,i) = MD_total;

      
%%%Rank Order Rearrangement


Ranked_AvgConnProf_acrossTime_Group1 = tiedrank( squeeze(AvgConnProf_acrossTime_Group1(b,:,i)) ); 
Ranked_AvgConnProf_acrossTime_Group2 = tiedrank( squeeze(AvgConnProf_acrossTime_Group2(b,:,i)) );  

  
  for j=1:Total_Target_ROIs
  
    MD_rank(j)=abs(Ranked_AvgConnProf_acrossTime_Group1(b,j,i)-Ranked_AvgConnProf_acrossTime_Group2(b,j,i));
  
  end 
  
    MD_rank_total=sum(MD_rank);
    Subject_RankOrderRearrangement(b,i) = MD_rank_total;


      
  end
  end

GroupAvg_Aggregate_Divergence = mean(Subject_Aggregate_Divergence,2);
GroupAvg_RankOrderRearrangement = mean(Subject_RankOrderRearrangement,2);        
GroupAvg_Entropy_Difference = mean(Subject_Entropy_Difference,2); 

%%%% Project the Subject-Level maps into an MNI space NIFTI file %%%%

imgRoot = '/path/to/seedmask';
mask = load_untouch_nii(sprintf('%s/Seed_Mask.nii.gz', imgRoot));

in_brain=find(mask.img(:)>0);
full_brain_space=mask.img.*0;
full_brain_space = double(full_brain_space);

for i=1:Subjects

%Aggregate Divergence Files
full_brain_space(in_brain)=Subject_Aggregate_Divergence(:,i) 
imgFile = sprintf('%s/Seed_Mask.nii.gz', imgRoot);
gm = load_untouch_nii(imgFile);
new_hdr=gm;
new_hdr.vol=full_brain_space;
new_hdr.img=full_brain_space;
flNm = num2str(i,'Voxelwise_Aggregate_Divergence_Subj_00%d');
save_untouch_nii(new_hdr,flNm);

%Rank Order Rearrangement Files
full_brain_space(in_brain)=Subject_RankOrderRearrangement(:,i) 
imgFile = sprintf('%s/Seed_Mask.nii.gz', imgRoot);
gm = load_untouch_nii(imgFile);
new_hdr=gm;
new_hdr.vol=full_brain_space;
new_hdr.img=full_brain_space;
flNm = num2str(i,'Voxelwise_RankOrder_Rearrangement_Subj_00%d');
save_untouch_nii(new_hdr,flNm);

%Entropy Shift Files
full_brain_space(in_brain)=Subject_Entropy_Difference(:,i) 
imgFile = sprintf('%s/Seed_Mask.nii.gz', imgRoot);
gm = load_untouch_nii(imgFile);
new_hdr=gm;
new_hdr.vol=full_brain_space;
new_hdr.img=full_brain_space;
flNm = num2str(i,'Voxelwise_Entropy_Shift_Subj_00%d');
save_untouch_nii(new_hdr,flNm);

end


  
  %Optional: CSVs containing value of connectivity profile metric for each voxel in the
  %mask
  writematrix(GroupAvg_Aggregate_Divergence, ['AggregateDivergence.csv'])
  writematrix(GroupAvg_Entropy_Difference, ['EntropyShift.csv'])
  writematrix(GroupAvg_RankOrderRearrangement, ['RankOrderRearrangement.csv'])
