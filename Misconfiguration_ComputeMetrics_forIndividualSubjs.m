%% Rather than compute connectivity profile metrics for one group versus
%% another group, compute connectivity profile metrics at the individual subject
%% level. Establish a set of normative average connectivity profiles, then compare
%% each subject in the clinical/treatment group to the normative average




clc;
clear;
close all;

softwareRoot = '/Applications/MATLAB_R2020a.app';;
%script requires download of NIFTI and ANALYZE image toolbox for Matlab
%(https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
addpath '/NIfTI_20140122'
addpath(genpath('/MATLAB'))'

%load CSV with subject IDs%
csvRoot = '...';
csvFile = sprintf('%s/SubjectIDs.csv', csvRoot);
fileId = fopen(csvFile);
csv = textscan(fileId, '%s', 'Delimiter', ',');
fclose(fileId);

%load CSV with seed IDs%
csv2Root = '...';
csv2File = sprintf('%s/Seeds.csv', csv2Root);
file2Id = fopen(csv2File);
csv2 = textscan(file2Id, '%s', 'Delimiter', ',');
fclose(file2Id);

%path to processed functional imaging data (z-scored, voxel-wise
%connectivity maps, one map for each ROI for each subject)
imgRoot = '.../ImagingData';

%load mask
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


intraGroup_1 = Y(1:33);   %Controls
intraGroup_2 = Y(34:end); %SUD

numVolIds_intraGroup_1 = length(intraGroup_1);
numVolIds_intraGroup_2 = length(intraGroup_2);

W1 = zeros(numVolIds_intraGroup_1,numVoxels);
Z1 = zeros(numSeeds,1);

W2 = zeros(numSeeds,numVoxels);
Z2 = zeros(numSeeds,1);

Z_final_MD = zeros(length(intraGroup_2),1);
Z_final_RankOrder = zeros(length(intraGroup_2),1);
Z_final_Entropy = zeros(length(intraGroup_2),1);


%create group average connectivity profiles for Control group

for j = 1:numSeeds
    seedId = csv2{:,1}{j,1};
    
    
W1 = zeros(numVolIds_intraGroup_1,numVoxels);

 for i = 1:numVolIds_intraGroup_1
    subjId = intraGroup_1{i,1};
   
    X1 = zeros(1, numVoxels);

    imgFile = sprintf('Controls/subj%s_r2zmap_%s.nii.gz', subjId, seedId);
   
    gm = load_nii(imgFile);
    gm = gm.img(maskIdx);
    X1(1, 1:numVoxels) = gm;

      
   W1(i,:) = X1;
    
  end
  
  W1 = mean(W1,2);
  
  Z1(j,:) = mean(W1, 1);
  Z1ranked = tiedrank(Z1);
  
  end

  %create individual connectivity profiles for each subject in SUD group,
  %compare each to Control group average
  
  for i = 1:numVolIds_intraGroup_2
    subjId = intraGroup_2{i,1};
   
    W2 = zeros(numSeeds,numVoxels);
    
    for j = 1:numSeeds
    seedId = csv2{:,1}{j,1};
    
    X2 = zeros(1, numVoxels);

    imgFile = sprintf('SUD/subj%s_r2zmap_%s.nii.gz', subjId, seedId);
   
    
    gm = load_nii(imgFile);
    gm = gm.img(maskIdx);
    X2(1, 1:numVoxels) = gm;

      
   W2(j,:) = X2; 
    
  end

   W2 = mean(W2,2);
  
  Z2ranked = tiedrank(W2);
  
  
  
  %compute the connectivity profile metrics at each voxel for the SUD subject
  for b=1:numVoxels
   
  %Entropy    
    h1= abs(Z1(:,b));
    p1 = bsxfun(@rdivide,h1,sum(h1,1));
    e1 = -nansum(p1.*log2(p1),1);
    Entropy1 = e1/log2(30);
      
    
    h2= abs(W2(:,b));
    p2 = bsxfun(@rdivide,h2,sum(h2,1));
    e2 = -nansum(p2.*log2(p2),1);
    Entropy2 = e2/log2(30);
    
  Z_final_Entropy(i,b) = abs(Entropy1-Entropy2);
  
  %Manhattan distance (aggregate divergence)
      MD1=abs(Z1(1,b)-W2(1,b));
      MD2=abs(Z1(2,b)-W2(2,b));
      MD3=abs(Z1(3,b)-W2(3,b));
      MD4=abs(Z1(4,b)-W2(4,b));
      MD5=abs(Z1(5,b)-W2(5,b));
      MD6=abs(Z1(6,b)-W2(6,b));
      MD7=abs(Z1(7,b)-W2(7,b));
      MD8=abs(Z1(8,b)-W2(8,b));
      MD9=abs(Z1(9,b)-W2(9,b));
      MD10=abs(Z1(10,b)-W2(10,b));
      MD11=abs(Z1(11,b)-W2(11,b));
      MD12=abs(Z1(12,b)-W2(12,b));
      MD13=abs(Z1(13,b)-W2(13,b));
      MD14=abs(Z1(14,b)-W2(14,b));
      MD15=abs(Z1(15,b)-W2(15,b));
      MD16=abs(Z1(16,b)-W2(16,b));
      MD17=abs(Z1(17,b)-W2(17,b));
      MD18=abs(Z1(18,b)-W2(18,b));
      MD19=abs(Z1(19,b)-W2(19,b));
      MD20=abs(Z1(20,b)-W2(20,b));
      MD21=abs(Z1(21,b)-W2(21,b));
      MD22=abs(Z1(22,b)-W2(22,b));
      MD23=abs(Z1(23,b)-W2(23,b));
      MD24=abs(Z1(24,b)-W2(24,b));
      MD25=abs(Z1(25,b)-W2(25,b));
      MD26=abs(Z1(26,b)-W2(26,b));
      MD27=abs(Z1(27,b)-W2(27,b));
      MD28=abs(Z1(28,b)-W2(28,b));
      MD29=abs(Z1(29,b)-W2(29,b));
      MD30=abs(Z1(30,b)-W2(30,b));
      
      MD_Total=MD1+MD2+MD3+MD4+MD5+MD6+MD7+MD8+MD9+MD10+MD11+MD12+MD13+MD14+MD15+MD16+MD17+MD18+MD19+MD20+MD21+MD22+MD23+MD24+MD25+MD26+MD27+MD28+MD29+MD30;
     
  
  
  Z_final_MD_sated(i,b) = MD_Total;
      
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
     
  
  
  Z_final_RankOrder(i,b) = MD_Total_Rank;
      
  end
  
  end
