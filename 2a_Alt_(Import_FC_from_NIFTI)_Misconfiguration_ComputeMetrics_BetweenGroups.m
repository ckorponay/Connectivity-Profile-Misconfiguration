%%Use if Script 1 was skipped and you want to import FC values directly from voxel-wise NIFTI files%%

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

%load CSV with subject IDs%
csvRoot = '...';
csvFile = sprintf('%s/SubjectIDs.csv', csvRoot);
fileId = fopen(csvFile);
csv = textscan(fileId, '%s', 'Delimiter', ',');
fclose(fileId);

%load CSV with ROI/seed IDs (i.e., the connectivity profile target regions)
csv2Root = '...';
csv2File = sprintf('%s/Seeds.csv', csv2Root);
file2Id = fopen(csv2File);
csv2 = textscan(file2Id, '%s', 'Delimiter', ',');
fclose(file2Id);

%path to processed functional imaging data (z-scored, voxel-wise
%connectivity maps, one map for each connectivity profile target ROI for each subject)
imgRoot = '.../ImagingData';

%load mask (mask region should encompass all voxels whose connectivity profiles will be evaluated)
mask = load_nii(sprintf('.../Mask.nii.gz'));
maskIdx = find(mask.img(:) > 0);

numVolIds = length(csv{:, 1});
numSeeds = length(csv2{:, 1});
numVoxels = length(maskIdx);

subjIds = cell(numVolIds, 1);
SeedsIds = cell(numSeeds, 1);



for i = 1:numVolIds
   
     Y(i, 1) = (csv{:, 1}(i));
    
end



P = 1

Z_final_MD = zeros(P,numVoxels);               %For storing the final, voxel-wise Aggregate Divergence values
Z_final_RankOrder = zeros(P,numVoxels);        %For storing the final, voxel-wise Rank Order Rearrangement values
Z_final_Entropy = zeros(P,numVoxels);          %For storing the final, voxel-wise Entropy Shift values


%%%% Set up each voxel's connectivity profile in each subject, and then create group-average connectivity profiles for each voxel


fprintf('\nStarting Iteration %d', p);

intraGroup_1 = Y(1:33);   %Controls            %e.g., for a control group of n=33
intraGroup_2 = Y(34:end); %SUD group    

numVolIds_intraGroup_1 = length(intraGroup_1);
numVolIds_intraGroup_2 = length(intraGroup_2);

W1 = zeros(numVolIds_intraGroup_1,numVoxels);
Z1 = zeros(numSeeds,numVoxels);

W2 = zeros(numVolIds_intraGroup_2,numVoxels);
Z2 = zeros(numSeeds,numVoxels);


%Group 1
for j = 1:numSeeds
    seedId = csv2{:,1}{j,1};

 for i = 1:numVolIds_intraGroup_1
    subjId = intraGroup_1{i,1};
   
    X1 = zeros(1, numVoxels);

    imgFile = sprintf('Controls/subj%s_r2zmap_%s.nii.gz', subjId, seedId);
   
    gm = load_nii(imgFile);
    gm = gm.img(maskIdx);
    X1(1, 1:numVoxels) = gm;

      
   W1(i,:) = X1;
    
 end
  
  %for getting the mean Z score of each ROI across all voxels
  Z1(j,:) = mean(W1, 1);
  %create rank ordered connectivity profile at each voxel
  Z1ranked = tiedrank(Z1);
  
end
  
  
  
%Group 2
  for j = 1:numSeeds
    seedId = csv2{:,1}{j,1};

 for i = 1:numVolIds_intraGroup_2
    subjId = intraGroup_2{i,1};
   
    X2 = zeros(1, numVoxels);

    imgFile = sprintf('SUD/subj%s_r2zmap_%s.nii.gz', subjId, seedId);
    
    gm = load_nii(imgFile);
    gm = gm.img(maskIdx);
    X2(1, 1:numVoxels) = gm;

      
   W2(i,:) = X2; 
    
 end
  
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
  
 
