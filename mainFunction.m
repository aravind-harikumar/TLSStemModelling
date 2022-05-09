function mainFunction()

% warning off; clearvars;
clear all; close all; rng(1);

% SET BASE FOLDER AND DEPENDENCEIES
baseFolder = '.'; % '.' corresponds to current folder
addpath(genpath(strcat(baseFolder,'/Matlab_depedencies/')));

% INPUT AND OUTPUT FOLDER PATHS
PARAM.InOutBasePath = strcat(baseFolder,'/LiDARDataSingleTrees/FGI_data/');
PARAM.InTLSDataFolder = 'plot_birch';
PARAM.OutTLSDataFolder = 'Cropped_TLS';
InTLSPlotNames = {'1062'}; % Input TLS plot data

%GLOBAL PARAMETERS AND PRINT FLAGS
PARAM.PerformTLSCropping = false;
PARAM.PerformStemParamterEstimation = true;
% Print flags
PARAM.PLOT_ON_GLOBAL_FLAG = false; % Shows Fig1: Entire Tree-level
PARAM.PLOT_ON_GLOBAL_FLAG_FIG2 = false; % Shows Fig2: Stem-Crossection
PARAM.PRINT_STEM = false; % Print stem on/off flag
PARAM.PRINT_TREE_CLOUD = false; % Print Raw LiDAR data flag
% Set Voxel Size
PARAM.VOXEL_SIZE = 0.01; % 0.01 = 1cm

% PERFORM TLS DATA EXTRACTION AND TREE PARAMETER ESTIMATION
for i=1:1:size(InTLSPlotNames,2)
    PARAM.In_data_name = InTLSPlotNames{i}; % TLS plot folder / file name
    disp(strcat('Currently working on:', {' '},num2str(PARAM.In_data_name)));
    % [PCData] = getInvertedCloud(inFilepath, OutFilePath, speciesFolder{i}, inParams);
    PCData = getStemPointsAndKnots(PARAM); % Detect 3D stem and branch-knots.
end

end
