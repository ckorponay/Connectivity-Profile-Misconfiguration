%Uses input from Script 1%

% Script for computing voxel-wise connectivity profile metrics (i.e., aggregate
% divergence, rank order misarrangement, entropy shift) between a clinical
% group and normal comparison group or between timepoints of the same group


clc;
clear;
close all;

softwareRoot = '/Applications/MATLAB_R2020a.app';;
%script requires download of NIFTI and ANALYZE image toolbox for Matlab
%(https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
addpath '/NIfTI_20140122'
addpath(genpath('/MATLAB'))


Final_Aggregate_Divergence = zeros(numVoxels);          %For storing the final, voxel-wise Aggregate Divergence values
Final_RankOrderRearrangement = zeros(numVoxels);        %For storing the final, voxel-wise Rank Order Rearrangement values
Final_Entropy_Difference = zeros(numVoxels);            %For storing the final, voxel-wise Entropy Shift values

GroupAvg_ConnProf_Group1=GroupAvg_ConnProf_Group1';
GroupAvg_ConnProf_Group2=GroupAvg_ConnProf_Group2';
  
%%%%%%%%%
%%Compute the 3 connectivity profile properties at each voxel, by comparing connectivity profiles across the group-averaged maps at each voxel

  for b=1:numVoxels
   
%%%Entropy Difference

    h1= abs(GroupAvg_ConnProf_Group1(:,b));
    p1 = bsxfun(@rdivide,h1,sum(h1,1));
    e1 = -nansum(p1.*log2(p1),1);
    Entropy_Group1 = e1/log2(Total_Target_ROIs);
      
    
    h2= abs(GroupAvg_ConnProf_Group2(:,b));
    p2 = bsxfun(@rdivide,h2,sum(h2,1));
    e2 = -nansum(p2.*log2(p2),1);
    Entropy_Group2 = e2/log2(Total_Target_ROIs);
    
    Final_Entropy_Difference(b) = abs(Entropy_Group1-Entropy_Group2);
  
%%%Aggregate Divergence (Manhattan Distance)

  for j=1:Total_Target_ROIs
  
    MD(j)=abs(GroupAvg_ConnProf_Group1(j,b)-GroupAvg_ConnProf_Group2(j,b));
  
  end 
  
    MD_total=sum(MD);
    Final_Aggregate_Divergence(b) = MD_total;

      
%%%Rank Order Rearrangement

  Ranked_GroupAvg_ConnProf_Group1 = tiedrank(GroupAvg_ConnProf_Group1);
  Ranked_GroupAvg_ConnProf_Group2 = tiedrank(GroupAvg_ConnProf_Group2);
  
  for j=1:Total_Target_ROIs
  
    MD_rank(j)=abs(Ranked_GroupAvg_ConnProf_Group1(j,b)-Ranked_GroupAvg_ConnProf_Group2(j,b));
  
  end 
  
    MD_rank_total=sum(MD_rank);
    Final_RankOrderRearrangement(b) = MD_rank_total;


      
  end
  
  %CSVs contain value of connectivity profile metric for each voxel in the
  %mask
  writematrix(Final_Aggregate_Divergence, ['AggregateDivergence.csv'])
  writematrix(Final_Entropy_Difference, ['EntropyShift.csv'])
  writematrix(Final_RankOrderRearrangement, ['RankOrderRearrangement.csv'])
