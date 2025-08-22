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


Z_final_MD = zeros(P,numVoxels);               %For storing the final, voxel-wise Aggregate Divergence values
Z_final_RankOrder = zeros(P,numVoxels);        %For storing the final, voxel-wise Rank Order Rearrangement values
Z_final_Entropy = zeros(P,numVoxels);          %For storing the final, voxel-wise Entropy Shift values



  
  %for getting the mean Z score of each ROI across all voxels
  Z2(j,:) = mean(W2, 1);
  %create rank ordered connectivity profile at each voxel
  Z2ranked = tiedrank(Z2);
  
  end
  
  
%%%%%%%%%
%%Compute the 3 connectivity profile properties at each voxel, by comparing connectivity profiles across the group-averaged maps at each voxel

%Set up is for connectivity profiles with 30 target regions (can be adjusted to any number)

  for b=1:numVoxels
   
  %Entropy    
    h1= abs(Z1(:,b));
    p1 = bsxfun(@rdivide,h1,sum(h1,1));
    e1 = -nansum(p1.*log2(p1),1);
    Entropy1 = e1/log2(30);
      
    
    h2= abs(Z2(:,b));
    p2 = bsxfun(@rdivide,h2,sum(h2,1));
    e2 = -nansum(p2.*log2(p2),1);
    Entropy2 = e2/log2(30);
    
  Z_final_Entropy(p,b) = abs(Entropy1-Entropy2);
  
  %Manhattan Distance (aggregate divergence)
      MD1=abs(Z1(1,b)-Z2(1,b));
      MD2=abs(Z1(2,b)-Z2(2,b));
      MD3=abs(Z1(3,b)-Z2(3,b));
      MD4=abs(Z1(4,b)-Z2(4,b));
      MD5=abs(Z1(5,b)-Z2(5,b));
      MD6=abs(Z1(6,b)-Z2(6,b));
      MD7=abs(Z1(7,b)-Z2(7,b));
      MD8=abs(Z1(8,b)-Z2(8,b));
      MD9=abs(Z1(9,b)-Z2(9,b));
      MD10=abs(Z1(10,b)-Z2(10,b));
      MD11=abs(Z1(11,b)-Z2(11,b));
      MD12=abs(Z1(12,b)-Z2(12,b));
      MD13=abs(Z1(13,b)-Z2(13,b));
      MD14=abs(Z1(14,b)-Z2(14,b));
      MD15=abs(Z1(15,b)-Z2(15,b));
      MD16=abs(Z1(16,b)-Z2(16,b));
      MD17=abs(Z1(17,b)-Z2(17,b));
      MD18=abs(Z1(18,b)-Z2(18,b));
      MD19=abs(Z1(19,b)-Z2(19,b));
      MD20=abs(Z1(20,b)-Z2(20,b));
      MD21=abs(Z1(21,b)-Z2(21,b));
      MD22=abs(Z1(22,b)-Z2(22,b));
      MD23=abs(Z1(23,b)-Z2(23,b));
      MD24=abs(Z1(24,b)-Z2(24,b));
      MD25=abs(Z1(25,b)-Z2(25,b));
      MD26=abs(Z1(26,b)-Z2(26,b));
      MD27=abs(Z1(27,b)-Z2(27,b));
      MD28=abs(Z1(28,b)-Z2(28,b));
      MD29=abs(Z1(29,b)-Z2(29,b));
      MD30=abs(Z1(30,b)-Z2(30,b));
      
      MD_Total=MD1+MD2+MD3+MD4+MD5+MD6+MD7+MD8+MD9+MD10+MD11+MD12+MD13+MD14+MD15+MD16+MD17+MD18+MD19+MD20+MD21+MD22+MD23+MD24+MD25+MD26+MD27+MD28+MD29+MD30;
     
  
  
  Z_final_MD(p,b) = MD_Total;
      
  %Rank Order Arrangement
      MD1=abs(Z1ranked(1,b)-Z2ranked(1,b));
      MD2=abs(Z1ranked(2,b)-Z2ranked(2,b));
      MD3=abs(Z1ranked(3,b)-Z2ranked(3,b));
      MD4=abs(Z1ranked(4,b)-Z2ranked(4,b));
      MD5=abs(Z1ranked(5,b)-Z2ranked(5,b));
      MD6=abs(Z1ranked(6,b)-Z2ranked(6,b));
      MD7=abs(Z1ranked(7,b)-Z2ranked(7,b));
      MD8=abs(Z1ranked(8,b)-Z2ranked(8,b));
      MD9=abs(Z1ranked(9,b)-Z2ranked(9,b));
      MD10=abs(Z1ranked(10,b)-Z2ranked(10,b));
      MD11=abs(Z1ranked(11,b)-Z2ranked(11,b));
      MD12=abs(Z1ranked(12,b)-Z2ranked(12,b));
      MD13=abs(Z1ranked(13,b)-Z2ranked(13,b));
      MD14=abs(Z1ranked(14,b)-Z2ranked(14,b));
      MD15=abs(Z1ranked(15,b)-Z2ranked(15,b));
      MD16=abs(Z1ranked(16,b)-Z2ranked(16,b));
      MD17=abs(Z1ranked(17,b)-Z2ranked(17,b));
      MD18=abs(Z1ranked(18,b)-Z2ranked(18,b));
      MD19=abs(Z1ranked(19,b)-Z2ranked(19,b));
      MD20=abs(Z1ranked(20,b)-Z2ranked(20,b));
      MD21=abs(Z1ranked(21,b)-Z2ranked(21,b));
      MD22=abs(Z1ranked(22,b)-Z2ranked(22,b));
      MD23=abs(Z1ranked(23,b)-Z2ranked(23,b));
      MD24=abs(Z1ranked(24,b)-Z2ranked(24,b));
      MD25=abs(Z1ranked(25,b)-Z2ranked(25,b));
      MD26=abs(Z1ranked(26,b)-Z2ranked(26,b));
      MD27=abs(Z1ranked(27,b)-Z2ranked(27,b));
      MD28=abs(Z1ranked(28,b)-Z2ranked(28,b));
      MD29=abs(Z1ranked(29,b)-Z2ranked(29,b));
      MD30=abs(Z1ranked(30,b)-Z2ranked(30,b));
      
      MD_Total_Rank=MD1+MD2+MD3+MD4+MD5+MD6+MD7+MD8+MD9+MD10+MD11+MD12+MD13+MD14+MD15+MD16+MD17+MD18+MD19+MD20+MD21+MD22+MD23+MD24+MD25+MD26+MD27+MD28+MD29+MD30;
     
  
  
  Z_final_RankOrder(p,b) = MD_Total_Rank;
      
  end
  
  %CSVs contain value of connectivity profile metric for each voxel in the
  %mask
  writematrix(Z_final_MD, ['MD.csv'])
  writematrix(Z_final_entropy, ['Entropy.csv'])
  writematrix(Z_final_RankOrder, ['RankOrder.csv'])
