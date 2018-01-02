clc
clear

Image_directory = {'/home/user03/wuxiang/FCMap/Front_Inf_L/'}; % Brain image directory
OutputName = '/home/user03/wuxiang/testcode/T_lIFG.nii'; % Output name and directory
SeedSeries = '/home/user03/wuxiang/Behavioral/score_inde_inter.txt'; % Behavioral data for correlation 
MaskFile = '/home/user03/wuxiang/Behavioral/GM_mask_fun.nii'; % Mask directory
ImageCell = {}; % Image covariates
cov = load('/home/user03/wuxiang/Behavioral/cov.txt'); 
TextCell = {cov}; % Text covariates
Config.Value = 1; % 1: One-sample T; 2: Two-sample T; 3: Paired T; 4: ANCOVA; 5: ANCOVA repeat; 6: Corr
Config.pValue = 0.01; % p voxel
Config.pCorrect = 0.05; % Corrected p
Config.correct = 2; % 1: FSL smoothest + Alphasim; 2: AFNI 3dClustSim 

stat_pipeline( Image_directory, OutputName, MaskFile, ImageCell, TextCell, Config);