%%%%% classification using SVM-RFE model %%%%
close all; clear all; clc 
global options_svm ; 
options_svm = '-q -s 1 -t 0 -b 1'; %%-set parameters for SVM 
load LTLE_RTLE_NC_148ROI_features.mat
%% read feature form stats
dir1 = pwd;
addpath([dir1 filesep 'function']);
FeatureName = {'NumVert';'SurfArea';'GrayVol';'ThickAvg';'ThickStd';'MeanCurv';'GausCurv';'FoldInd';'CurvInd'};
group1 = 'RTLE';
group2 = 'NC';
Feature_ID = 0;
if Feature_ID ~= 0
    CurrFeatureName = FeatureName{Feature_ID};
else
    CurrFeatureName = 'All features';
end
flagNorm = 1; % normalize the feature vetor
options.c = 2;  % Class count
options.ITER_NUM = 10;
options.r = 1.0/2;  %% step size or control factor;
options.Thres = 0 ;  %% threshold for selectiing features
options.options_svm = ['-q -s 1 -t 0 -b 1']; 
options.CrossValidMode = 1; 
options.CrossValidPara = 0.10;

%% 2. load feature files
% load([group1 'features.mat']);
% load([group2 'features.mat']);
Feature_C1 = eval([group1 'features']);
Feature_C2 = eval([group2 'features']);
if Feature_ID ~= 0
    Feature_C1 = Feature_C1(:,:,Feature_ID);
    Feature_C2 = Feature_C2(:,:,Feature_ID); 
    if flagNorm == 1    % normalize feature
         Feature_C1 = zscore(Feature_C1,0,2);
         Feature_C2 = zscore(Feature_C2,0,2);
    end
else   % Select all features
    Feature_C1_All = [];
    Feature_C2_All = [];
    [SubjectNum,BrainRegionNum,FeatureNum] = size(Feature_C1);
    Featureindex1 = [2 3 4 6];
    for i =1:4
%     for i=1: FeatureNum
%         CurrFeature_C1 = Feature_C1(:,:,i);
%         CurrFeature_C2 = Feature_C2(:,:,i);
        CurrFeature_C1 = Feature_C1(:,:, Featureindex1(i));  %add by Lai 2015/12/22
        CurrFeature_C2 = Feature_C2(:,:, Featureindex1(i)); 
        if flagNorm == 1
           CurrFeature_C1 = zscore(CurrFeature_C1,0,2);
           CurrFeature_C2 = zscore(CurrFeature_C2,0,2);
        end
        Feature_C1_All = [Feature_C1_All CurrFeature_C1];
        Feature_C2_All = [Feature_C2_All CurrFeature_C2];
    end
    Feature_C1 = Feature_C1_All;
    Feature_C2 = Feature_C2_All;
end

Feature = [Feature_C1;Feature_C2];
PosLabel = 1; 
NegLabel = -1;
Label = ones(size(Feature_C1,1)+size(Feature_C2,1),1); 
Label(1:size(Feature_C1,1)) = PosLabel;
Label(size(Feature_C1,1)+1:end)=NegLabel;
% Feature=mapminmax(Feature,0,1);
Feature=zscore(Feature);

%% Feature selection -- main function --
[FeatureRank,PerformanceMatrix,ScoreMatrix,Num_select,FeatureIndexSortBest,dec,Label_roc] = FS_Rank_freesurfer(Feature,Label);
FeatureRank = fliplr(FeatureRank);  % sort the rank form important to unimportant
% FeatureRank is the featureIndex in different brain region,from useless to useful.
% PerformanceMatrix is the performance in cross validation from much dimension to less

% load('Brainregion_aparc.mat'); 
% region_Rank = region(FeatureRank); % region from important to unimportant
ACCArray = PerformanceMatrix(1,:); 
TPRArray = PerformanceMatrix(2,:); 
TNRArray = PerformanceMatrix(3,:); 

savename = [group1 '-' group2 '-SVM-RFE-' CurrFeatureName];
save(savename);
plotm
plotROC;