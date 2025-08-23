%Uses input from Script 2%

% Script for computing voxel-wise connectivity profile metrics (i.e., aggregate
% divergence, rank order misarrangement, entropy shift) between different timepoints of randomly
% sampled subsets of a large normative dataset


clc;
clear;
close all;

softwareRoot = '/Applications/MATLAB_R2020a.app';;
%script requires download of NIFTI and ANALYZE image toolbox for Matlab
%(https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
addpath '/NIfTI_20140122'
addpath(genpath('/MATLAB'))

Num_Permutations=10000

Normative_Aggregate_Divergence = zeros(Num_Permutations,Seed_Voxels);          %For storing the final, voxel-wise Aggregate Divergence values
Normative_RankOrderRearrangement = zeros(Num_Permutations,Seed_Voxels);        %For storing the final, voxel-wise Rank Order Rearrangement values
Normative_Entropy_Difference = zeros(Num_Permutations,Seed_Voxels);            %For storing the final, voxel-wise Entropy Shift values

for P=1:Num_Permutations


  
%%%%%%%%%
%%Compute the change between time points in each of the 3 connectivity profile properties at each voxel in each individual 

for i=1:Subjects

  for b=1:Seed_Voxels
   
%%%Entropy Difference

    h1= abs(AvgConnProf_acrossTime_NormativeData(b,:,i));
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

  Ranked_AvgConnProf_acrossTime_Group1 = tiedrank(AvgConnProf_acrossTime_Group1);
  Ranked_AvgConnProf_acrossTime_Group2 = tiedrank(AvgConnProf_acrossTime_Group2);
  
  for j=1:Total_Target_ROIs
  
    MD_rank(j)=abs(Ranked_AvgConnProf_acrossTime_Group1(b,j,i)-Ranked_AvgConnProf_acrossTime_Group1(b,j,i));
  
  end 
  
    MD_rank_total=sum(MD_rank);
    Subject_RankOrderRearrangement(b,i) = MD_rank_total;


      
  end
  end

GroupAvg_Aggregate_Divergence = mean(Subject_Aggregate_Divergence,2);
GroupAvg_RankOrderRearrangement = mean(Subject_RankOrderRearrangement,2);        
GroupAvg_Entropy_Difference = mean(Subject_Entropy_Difference,2); 

  
  %CSVs contain value of connectivity profile metric for each voxel in the
  %mask
  writematrix(GroupAvg_Aggregate_Divergence, ['AggregateDivergence.csv'])
  writematrix(GroupAvg_Entropy_Difference, ['EntropyShift.csv'])
  writematrix(GroupAvg_RankOrderRearrangement, ['RankOrderRearrangement.csv'])
