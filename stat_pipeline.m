function [ CorrectImage ] = stat_pipeline( S, OutputName, MaskFile, ImageCell, TextCell, Config, SeedSeries)
%   This function is voxel-wised statistical pipeline including statistical and correction process.
%   Before using this function, make sure spm, REST, DPABI, FSL or AFNI have already added to the path.
% ------------------------------------------------------------------------------------------------------
%   Input:
%       S          - Image directory
%       OutputName - Output name and directory
%       MaskFile   - Mask directory
%       ImageCell  - Image covariates
%       TextCell   - Text covariates
%       Config     - Configuration for p values and corrrection method
%       SeedSeries - Behavioral data
%   Output:
%       Corrected_xxxxx.nii Cluster image after correction
%   By Ma Junji 20171024
%   Version = 1.0
%------------------------------------------------------------------------------------------------------

% Read correlation text or not
if nargin < 7
    SeedSeries = '';
else
    BehavioralData = load(SeedSeries);
end

% Define statictical method
switch Config.Value
    case 1 %One-Sample
        y_TTest1_Image(S, OutputName, MaskFile, ImageCell, TextCell, 0);
        testFlag = 'T';
    case 2 %Two-Sample
        y_TTest2_Image(S, OutputName, MaskFile, ImageCell, TextCell);
        testFlag = 'T';
    case 3 %Paired
        y_TTestPaired_Image(S, OutputName, MaskFile, ImageCell, TextCell);
        testFlag = 'T';
    case 4 %ANCOVA
        y_ANCOVA1_Image(S, OutputName, MaskFile, ImageCell, TextCell);
        testFlag = 'F';
        %         MC_list = {'None';'tukey-kramer';'lsd';'bonferroni';'dunn-sidak';'scheffe';}; %YAN Chao-Gan, 151127. Add multiple comparison test for ANCOVA
        %         MC_Value = get(handles.popupmenuMC,'Value');
        %         MC_type = MC_list{MC_Value};
        %         y_ANCOVA1_Multcompare_Image(S, OutputName, MaskFile, ImageCell, TextCell, MC_type);
    case 5 %ANCOVA Repeat
        y_ANCOVA1_Repeated_Image(S, OutputName, MaskFile, ImageCell, TextCell);
        testFlag = 'F';
    case 6 %Corr
        y_Correlation_Image(S, BehavioralData, OutputName, MaskFile, ImageCell, TextCell);
        testFlag = 'R';
end

% Change to the output directory and load stat map
[Path, Name, Ext]=fileparts(OutputName);
cd(Path);
ResultHeader = spm_vol([Name, Ext]);
[ResultImage, xxx] = spm_read_vols(ResultHeader);
DOF = CheckDf(ResultHeader);
pthr = Config.pValue;
athr = Config.pCorrect;

switch Config.correct
    case 1
        % Estimate FWHM using FSL
        Residual_Name = [Name,'_Residual','.nii'];
        if testFlag ~= 'F'
            smoothest_cmd = ['smoothest -m ',MaskFile,' -r ',Residual_Name,' -d ',num2str(DOF),' -V|grep "FWHMx" > FWHM.txt'];
        else
            smoothest_cmd = ['smoothest -m ',MaskFile,' -r ',Residual_Name,' -d ',num2str(DOF(2)),' -V|grep "FWHMx" > FWHM.txt'];
        end
        unix(smoothest_cmd);
        fid = fopen('FWHM.txt');
        FWHM_cell = textscan(fid,'%s');
        fclose(fid);
        FWHM = [str2double(FWHM_cell{1}{15}),str2double(FWHM_cell{1}{19}),str2double(FWHM_cell{1}{23})];
        % Calculate cluster size (Alphasim)
        rmm = 5;
        iter = 1000;
        outname = ['Alphasim_',Name];
        rest_AlphaSim(MaskFile, Path, outname, rmm, FWHM, pthr, iter);
        [~,~,~,~,~,alpha] = textread([outname,'.txt'],'%f %f %f %f %f %f','headerlines',22);
        Index = find(alpha < athr);
        ClusterSize = Index(1);
        
    case 2
        % Estimate acf using AFNI
        Residual_Name = [Name,'_Residual','.nii'];
        FWHM_cmd = ['3dFWHMx -mask ',MaskFile,' -dset ',Residual_Name,' -detrend -geom -acf > ACF.txt'];
        unix(FWHM_cmd);
        % Calculate cluster size
        ACF_matrix = load('ACF.txt');
        ACF_a = ACF_matrix(2,1); ACF_b = ACF_matrix(2,2); ACF_c = ACF_matrix(2,3);
        [Img_x, Img_y, Img_z] = size(ResultImage);
        ImgSize = ['-nxyz ',num2str(Img_x),' ',num2str(Img_y),' ',num2str(Img_z),' '];
        Voxel_x = -ResultHeader.mat(1,1); Voxel_y = ResultHeader.mat(2,2); Voxel_z = ResultHeader.mat(3,3);
        VoxelSize = ['-dxyz ',num2str(Voxel_x),' ',num2str(Voxel_y),' ',num2str(Voxel_z),' '];
        AcfValue = ['-acf ',num2str(ACF_a),' ',num2str(ACF_b),' ',num2str(ACF_c),' '];
        ClustSim_cmd = ['3dClustSim ',ImgSize,VoxelSize,AcfValue, '-mask ',MaskFile,' -pthr ',num2str(pthr),' -athr ',num2str(athr),' -NN 123 -iter 10000 -prefix ',Name];
        unix(ClustSim_cmd);
        [~,ClusterSize] = textread([Name,'.NN2_bisided.1D'],'%f %f','headerlines',8);
end

% Set significant threshold
if (testFlag == 'T')
    Thrd = tinv(1 - pthr/2,DOF); % Two-tailed
    % Thrd=tinv(1 - pthr,DOF); % One-tailed
elseif (testFlag == 'F')
    Thrd =finv(1-pthr,AConfig.Df.Ftest(1),AConfig.Df.Ftest(2));
elseif (testFlag == 'R')
    TRvalue=tinv(1 - pthr/2,DOF); % Two-tailed
    Thrd=sqrt(TRvalue^2/(DOF+TRvalue^2));
end
SigIndex = find(ResultImage>=Thrd|ResultImage<=-Thrd);
[I J K] = ind2sub(size(ResultImage),SigIndex);
ThrdImage = zeros(size(ResultImage));
ThrdImage(SigIndex) = ResultImage(SigIndex);

% Set cluster size and save clusters
Overlay.ClusterSizeThrd = ClusterSize;
Overlay.ClusterConnectivityCriterion = 18;
ThrdImage = ThrdOverlayCluster(Overlay,ThrdImage);
ResultHeader.dt(1) = 16;
OutHeader = ResultHeader;
cd(Path);
OutHeader.fname = ['corrected_',Name,'.nii'];
OutHeader = spm_write_vol(OutHeader,ThrdImage);

CorrectImage = ThrdImage;
end

function DOF=CheckDf(Header) % get degree of freedom
DOF = 0;
if isfield(Header,'descrip')
    headinfo = Header.descrip;
    if ~isempty(strfind(headinfo,'{T_['))
        Tstart=strfind(headinfo,'{T_[')+length('{T_[');
        Tend=strfind(headinfo,']}')-1;
        testDf = str2num(headinfo(Tstart:Tend));
    elseif ~isempty(strfind(headinfo,'{F_['))
        Tstart=strfind(headinfo,'{F_[')+length('{F_[');
        Tend=strfind(headinfo,']}')-1;
        testDf = str2num(headinfo(Tstart:Tend));
    elseif ~isempty(strfind(headinfo,'{R_['))
        Tstart=strfind(headinfo,'{R_[')+length('{R_[');
        Tend=strfind(headinfo,']}')-1;
        testDf = str2num(headinfo(Tstart:Tend));
    elseif ~isempty(strfind(headinfo,'{Z_['))
        testFlag='Z';
        Tstart=strfind(headinfo,'{Z_[')+length('{Z_[');
        Tend=strfind(headinfo,']}')-1;
        testDf = str2num(headinfo(Tstart:Tend));
    end
    DOF = testDf;
end
end

function Result =ThrdOverlayCluster(Overlay, AVolume)
%Threshold for Cluster Size or calculated cluster size from cluster-Raidus %This function must be called after thresholding the Value already!
%"AConfig.Overlay.VolumeThrd" must be thresholded before this function is called!
Result =AVolume;
%Raidus has Priority
% 	if AConfig.Overlay.ClusterRadiusThrd~=0, %Raidus(mm)
% 		%Calcute the cluster size according to the Raidus(mm)
% 		AConfig.Overlay.Header.Origin=AConfig.Overlay.Origin;
%         maskROI =rest_SphereROI( 'BallDefinition2Mask' , sprintf('ROI Center(mm)=(0, 0, 0); Radius=%g mm.', AConfig.Overlay.ClusterRadiusThrd), size(AConfig.Overlay.Volume), AConfig.Overlay.VoxelSize, AConfig.Overlay.Header);
% 		AConfig.Overlay.ClusterSizeThrd =length(find(maskROI));
% 		AConfig.Overlay.ClusterRadiusThrd =0;
% 	end
if Overlay.ClusterSizeThrd==0, return; end
[theObjMask, theObjNum]=bwlabeln(Result,Overlay.ClusterConnectivityCriterion); 	%[theObjMask, theObjNum]=bwlabeln(Result);
for x=1:theObjNum,
    theCurrentCluster = theObjMask==x;
    if length(find(theCurrentCluster))<Overlay.ClusterSizeThrd,
        Result(logical(theCurrentCluster))=0;
    end
end
end
