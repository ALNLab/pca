function [permPval,bootPval, newCompOut] = componentSig(inData,varargin)
%% [permPval,~] = componentSig(behavioralData,outputDir);
%
%% Read-me:
% This function will run pca or probabilistic pca on the matrix in inData
% (observations are rows, variables are columns)and generate some pvalues
% for each resulting component, as well as each loading within each
% component. 
% 
% Permutations of the data are taken to assess components significance
% (using only singular values from SVD, not full PCA coefficients), and
% bootstrap analysis is taken to assess loading significance (by looking at
% loading magnitude). Each permutations/bootstrap result is saved in the
% directory supplied in outDir. To facilitate longer analysis, all
% permutation/bootstrap results are split into bins that are saved as
% unique .mat files appended by their bin number. Each bin contains a
% matFiles variable that contains all .mat filepaths for the analysis that
% was specified. Significance testing with these permutation/bootstrap
% results is performed with the nonparamSig function which reads in
% matFiles and concatenates your data for you. It has some other handy
% features like generating z-scores (see nonparamSig help).
%
% You can choose to run only the permutation analysis by setting the number
% of permutations to be performed in numPerm to some number greater than
% zero and the number of bootstrap to be taken in numBoot to zero. The
% reverse can be done to only run bootstrap analysis without permutation.
%
% The output format of the script makes it such that you can do one
% analysis, then return to the other while reading in existing data. This
% is also useful if you are doing a longer analysis and you quit the
% script. If you rerun the same command, the script will pick up where it
% left off. When switching between analyses you can change
% overrideVars.Switch to overwrite existing data in outDir/results. For
% example, if you want to save both analyses at the end or if you want to
% overwrite the parameters saved in your results for a different analysis
% (e.g. switching from bootstrap to permutation).
%
% Either pca or probabilistic pca (ppca) will be performed depending on
% whether pcaType is set to 'normal' or 'probabilistic'. See documentation
% for matlab's ppca function for more details. If you have missing data,
% ppca is recommended and the script will warn you of this. The
% ppca_with_svd function 
% is used for permutation analysis with ppca. It is a version of matlab's
% ppca function that saves an SVD analysis immediately after data
% interpolation. If pca is performed the number of components will be equal
% to the number of variables (i.e. cols in inData) minus 1. If ppca is
% performed it will most likely be the number of variables OR number of
% subjects depending on which is smaller (minus 1). Ppca may fail to
% converge on this number of components so an iterative approach is
% implemented (i.e. you may have even fewer components as a result). For
% bootstrap analysis you can set a fixed number of components in the
% variable bootNumComp. 
%
% Bootstrap analysis can be run over scores as well as loadings. If so, you
% need to change bootVar to 'scores'. bootVar can be set to 'loadings' as
% well, or alternatively set it to 'both' to use the bootstrap appraoch to
% get p-values for both loadings and scores.  The script can also check
% for normality of bootstrapped distributions (this is useful to ensure
% that the values you get are proper z-scores (see below). 
%
% By default, both permutation and bootstrap analysis will center your
% data, but centerSwitch may be set to 'false' in order to use raw data for
% ONLY the permutation analysis. 
%
% Some added features : if projectedData.Switch is set to 'true', then
% projectedData.data can be a second matrix of the same size that will be
% projected into the first matrice's space using the loadings of the first
% matrix, thereby giving you scores of the second matrix in the space of the
% first matrix. This will be done through projectPCA script. That standalone
% script is also capable of generating some biplots of both scores but there
% is no option to do this in the current script.
%
%% Required Inputs:
% inData : [observations x variables] matrix
%
%% Optional Inputs:
% These must be added in order following inData:
%
% -outDir: directory in which subdirectories for each analysis will be
% placed
%
% -numPerm : number of permutations to run. If set to 0, then permutation
% analysis will be skipped.
%
% -rotateSwitch : If set to 'true', this will use orthogonal procrustes to
% rotate your components. This can also be done manually through nonparamSig
% but is only feasible if you save out loadings with your analysis. In the
% case of permutation analysis a manual procrustes script is applied (based
% on the script used in McIntosh's fMRI Partial Least Squares suite). In the
% case of the bootstrap analysis internal matlab scripts are used to rotate
% loadings. If set to 'false' no rotation is performed.
%
% - numBoot : number of bootstrap samples to take. If set to 0, then
% bootstrap analysis will be skipped.
%
% - bootVar : if set to 'loadings' bootstrap analysis will be run over
% loadings. if set to 'scores' bootstrap analysis will run over scores. If
% this is set to 'both' then the analysis will be performed for both scores
% and loadings. In this case, the first column in permPval will correspond
% to significance of loadings and the second column will correspond to
% scores.
%
% - bootNumComp : If this variable is empty then maximum number of
% components will
% be used for analysis. Recommended you have some a priori number of
% components otherwise this will be very time consuming (i.e. run
% permutatins first). In that case, just set this to equal some number.
%
% - normalCheckSwitch : If set to 'true', then this will assess each
% bootstrapped
% distribution for normality. This will significantly slow down the analysis
% but it will ensure the p-values are accurate (i.e. that we are generating
% z-scores; see methodology below). Otherwise, set to 'false'.
%
% - projectedData.Switch : If set to 'true', then this will project a
% second data matrix stored in projectedData.data into the pca space of
% inData using the loadings of inData (i.e. you will get scores). If
% bootstrap is greater than zero, then the script will assess projected
% score significance. Otherwise, set this to 'false'.
%
% - parallelSwitch : If set to 'true' will run all permutations/bootstraps
% within a bin in parallel. This may crash your computer if inData is too
% large. Otherwise, set to 'false'.
%
% - overrideVars.Switch : If set to 'all' this will overwrite the variables
% found in already generated data (in your outDir). If overrideVars.Switch
% is set to 'some' this will overwrite the variables found in already
% generated data that was found by the script, but this will only be done
% for variables supplied in overrideVars.vars. These variables must be in a
% cell where each row is a seperate variable to overwrite. For example,
% overrideVars.vars = {'numBoot';'numPerm'};. Otherwise, set
% overrideVars.Switch to 'false'.
%
% - centerSwitch : If set to 'true' this will center the data before
% applying SVD for permutation analysis. Otherwise, set to 'false'. NOTE:
% bootstrap analysis will still center the data and use the covariance
% matrix of inData to generate pca results.
%
%% Default Parameters:
% outDir : working directory
% numPerm : 500
% rotateSwitch : 'false'
% numBoot : 100
% bootVar : 'both'
% bootNumComp : ''
% normalCheckSwitch : 'false'
% projectedData.Switch : 'false'
% parallelSwitch : 'false'
% overrideVars.Switch : 'false'
% centerSwitch : 'true'
% pcaType : 'normal'
% rotateReal : 'false';
%
%% Example of call:
% using all input arguments: 
% [permPval,bootPval] = componentSig(inData,outDir,numPerm,rotateSwitch,numBoot,bootVar,bootNumComp,normalCheckSwitch,projectedData.Switch,parallelSwitch,overrideVars.Switch,centerSwitch,pcaType,rotateReal)
% run script in default:
% [permPval,bootPval] = componentSig(inData)
% run 10 permutations and rotate the SVD values according to the real data:
% [permPval,bootPval] = componentSig(inData,'/Users/ateghipc/OneDrive - personalmicrosoftsoftware.uci.edu/PT/pca/connectivity2',10,'true')
% run no permutations, but run 100 bootstrap samples.
% [permPval,bootPval] = componentSig(inData,'/Users/ateghipc/OneDrive - personalmicrosoftsoftware.uci.edu/PT/pca/connectivity2',0,'true',100,'loadings',3,'false','false','false','false','true')
%
%% Outputs:
% permPval : nx1 vector where rows correspond to components. Values are
% p-Values and represent significance of components based on permutation
% analysis.
% bootPval : rows correspond to loadings and columns to components. Values
% are p-values and represent significance of loadings based on bootstrap
% analysis.
%
%% Anatomy of stored mat files:
% realData.mat
% This contains several variables. First, all of the setup arguments will be
% saved for reference. If you asked to run permutations, this will contain
% svd output. Singular values are in S. The diagonal of singular values are
% in SD. The left hand singular values are in U and the right hand singular
% values are in V. If you asked to run bootstrap, this will also contain pca
% output. scores corresponds to eigenvectors, coeff are loadings, mu is the
% mean of each column of inData for centering. For other variables see pca
% documentation in matlab. If you asked to project data, projectedScores
% will correspond to your new new matrix in pca space.
%
% Perms_BinSize_BinNumber.mat
% matFiles is just a cell vector where each row corresponds to the full path
% of each of these mat files. The full size is determined by the full number
% of bins the script calculated it should run. It will be complete up to the
% current BinNumber. This variable facilitates calling of nonparamSig. It's
% saved out for each BinNumber in case script crashes. Then you can just
% feed this into nonparamSig outside of the current script.
% Sigma contains the permuted singular values. Rows correspond to components
% and columns to the permutation.
%
% Bootstrap_BinSize_BinNumber.mat
% See above for matFiles explanation. dataBoot columns correspond to
% components. The first X number of columns corresponds to the number of
% components from the first bootstrapped sample (corresponding to loadings
% or scores (depending on which of the two you made the script calculate).
%
%% Methodology
% Based on Mcintosh & Lobaugh (2004) published in Neuroimage.
% permutations : randomization of rows in inData. Component significance is
% assesed by evaluating permuted singular values.
% bootstrap : bootstrapped rows in inData.
%
%% Alex Notes:
% - PARALLEL BOOTSTRAP DOES NOT WORK
% - Make SVD procrustes a seperate script
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Alex Teghipco (alex.teghipco@uci.edu) Last update 12/1/17

%defaults
numBoot = 0; %taken from Mcintosh & Lobaugh (2004)
numPerm = 500; %taken from Mcintosh & Lobaugh (2004)
outDir = pwd;
bootNumComp = '';
rotateSwitch = 'true'; %taken from Mcintosh & Lobaugh (2004)
normalCheckSwitch = 'false';
projectedData.Switch = 'false';
parallelSwitch = 'false';
overrideVars.Switch = 'false';
bootVar = 'both'; %coeff from matlab pca is loadings
centerSwitch = 'true';
pcaType = 'normal';
rotateReal = 'false';

%Check if you asked for something other than default settings
fixedargn = 1;
if nargin > (fixedargn + 0)
    if ~isempty(varargin{1})
        outDir = varargin{1};
    end
end
if nargin > (fixedargn + 1)
    if ~isempty(varargin{2})
        numPerm = varargin{2};
    end
end
if nargin > (fixedargn + 2)
    if ~isempty(varargin{3})
        rotateSwitch = varargin{3};
    end
end
if nargin > (fixedargn + 3)
    if ~isempty(varargin{4})
        numBoot = varargin{4};
    end
end
if nargin > (fixedargn + 4)
    if ~isempty(varargin{5})
        bootVar = varargin{5};
    end
end
if nargin > (fixedargn + 5)
    if ~isempty(varargin{6})
        bootNumComp = varargin{6};
    end
end
if nargin > (fixedargn + 6)
    if ~isempty(varargin{7})
        normalCheckSwitch = varargin{7};
    end
end
if nargin > (fixedargn + 7)
    if ~isempty(varargin{8})
        projectedData.Switch = varargin{8};
    end
end
if nargin > (fixedargn + 8)
    if ~isempty(varargin{9})
        parallelSwitch = varargin{9};
    end
end
if nargin > (fixedargn + 9)
    if ~isempty(varargin{10})
        overrideVars.Switch = varargin{10};
    end
end
if nargin > (fixedargn + 10)
    if ~isempty(varargin{11})
        centerSwitch = varargin{11};
    end
end
if nargin > (fixedargn + 11)
    if ~isempty(varargin{12})
        pcaType = varargin{12};
    end
end
if nargin > (fixedargn + 12)
    if ~isempty(varargin{13})
        rotateReal = varargin{13};
    end
end

if sum(sum(isnan(inData))) > 0 && strcmp(pcaType,'normal') == 1
    warning('You have some nan values in your dataset ... I recommend running probabilistic pca for interpolation of missing data (see documentation)')
end

testLicense = license('test', 'Statistics_toolbox');
if testLicense == 0
    warn = warndlg('You do not have the necessary toolbox to run bootstrap analysis!');
end

%% Startup
h = waitbar(0,'Starting up ...');
outReal = [outDir '/Results/realData.mat'];
%Check if you've already run the analysis previously (real data)...
if exist(outReal,'file') == 2
    disp('Found an existing file with results ... loading it up now ...')
    switch overrideVars.Switch
        case 'all'
            load(outReal,'-regexp','^(?!numPerm$|bootNumComp$|numBoot$|outDir$|rotateSwitch$|normalCheckSwitch$|parallelSwitch$|bootVar$|h$).');
        case 'some'
            for i = 1:size(overrideVars.vars,1)+1
                if i == 1
                    overrideVars.vars{i} = ['^(!' overrideVars.vars{i} '$|'];
                end
                if i == size(overrideVars.vars,1)+1
                    overrideVars.vars{i} = 'h$).';
                end
                if i ~= 1 && i ~= size(overrideVars.vars,1)
                    overrideVars.vars{i} = [overrideVars.vars{i} '$|'];
                end
            end
            overVars = strjoin(overrideVars.vars,'');
            load(outReal,'-regexp',overVars);
        case 'false'
            load(outReal)
    end
    %h = waitbar(0,'Working on real data ...');
else
    %If you haven't, run pca. Seperately center the data and run SVD
    %(but only if you want to run the permutation analysis).
    mkdir([outDir '/Results'])
    %waitbar(0,h,'Working on real data ...');
end

waitbar(0,h,'Working on real data ...');


if exist('coeff','var') == 0 && numBoot > 0
    disp('Running pca on real data... ')
    if isempty(bootNumComp) == 1 && numBoot > 0
        disp('Bootstrap will use all available components ...')
        disp('Running pca ...')
        switch pcaType
            case 'normal'
                [coeff,score,latent,tsquared,explained,mu] = pca(inData);
                bootNumComp = size(score,2);
            case 'probabilistic'
                warning('Using maximum number of components ideally possible ... this may not work with ppca')
                bootNumComp = min(size(inData));
                tic
                    while exist('coeff','var') == 0
                        try
                            [~, ~, ~, coeff, score, pcvar, mu, v, ~] = ppca_with_svd(inData,bootNumComp);
                        catch
                            probAttempt = probAttempt - 1;
                            el = toc;
                            mins = 10;
                            if el > mins*60
                                error('There is an error with your data such that the ppca_with_svd function cannot be executed')
                            end
                        end
                    end
        end
    else
        if isempty(bootNumComp) == 0 && numBoot > 0
            switch pcaType
                case 'normal'
                [coeff,score,latent,tsquared,explained,mu] = pca(inData,'NumComponents',bootNumComp);
                case 'probabilistic'
                [coeff,score,pcvar,~,~,~] = ppca(inData,bootNumComp);
            end
        end
    end
    if exist('explained','var') == 0
        explained =100*pcvar/sum(pcvar);
    end
    switch rotateReal
        case 'false'
            disp('Not rotating real data ... ')
        case 'varimax'
            disp('Performing varimax rotation ... ')
            save([outDir '/Results/nonRotatedPCA.mat'],'coeff','score','explained')
            [coeff,rotateMat] = rotatefactors(coeff,'Method','varimax');
            score = score*rotateMat;
        case 'parsimax'
            disp('Performing parsimax rotation ... ')
            save([outDir '/Results/nonRotatedPCA.mat'],'coeff','score','explained')
            [coeff,rotateMat] = rotatefactors(coeff,'Method','parsimax');
            score = score*rotateMat;
        case 'promax'
            disp('Performing promax rotation ... ')
            save([outDir '/Results/nonRotatedPCA.mat'],'coeff','score','explained')
            [coeff,rotateMat] = rotatefactors(coeff,'Method','promax');
            score = score*rotateMat;
    end
end
if exist('S','var') == 0 && numPerm > 0
    disp('Running SVD on real data ...')
    switch centerSwitch
        case 'true'
            switch pcaType
                case 'normal'
                    inDataCentered = bsxfun(@minus,inData,mean(inData));
                    [U, S, V] = svd(inDataCentered,'econ');
                case 'probabilistic'
                    inDataCentered = bsxfun(@minus,inData,mean(inData,'omitnan'));
                    probAttempt = min(size(inDataCentered));
                    tic
                    while exist('U','var') == 0
                        try
                            [U, S, V, ~, ~, ~, ~, ~, ~] = ppca_with_svd(inDataCentered,probAttempt);
                        catch
                            probAttempt = probAttempt - 1;
                            el = toc;
                            mins = 10;
                            if el > mins*60
                                error('There is an error with your data such that the ppca_with_svd function cannot be executed')
                            end
                        end
                    end
            end
        case 'false'
            switch pcaType
                case 'normal'
                    [U, S, V] = svd(inData,'econ');
                case 'probabilistic'
                    probAttempt = min(size(inData));
                    while exist('U','var') == 0
                        try
                            [U, S, V, ~, ~, ~, ~, ~, ~] = ppca_with_svd(inData,probAttempt);
                        catch
                            probAttempt = probAttempt - 1;
                        end
                    end
            end
            
    end
    SD = diag(S);
    permNumComp = size(SD,1);
end

%If you asked to project some data, do that now as well...
%if exist('projectedData','var') == 0
    switch projectedData.Switch
        case 'true'
            disp('Projecting second matrix into first matrix pca space ... ')
            projectedScore = projectPCA(inData,'',projectedData.data,coeff,'false'); %projected score is ~V
            origData = inData;
            disp('Treating projected data matrix as primary input matrix ... ')
            inData = projectedData.data;
            savex(outReal,'inData','origData')
        case 'false'
            savex(outReal,'inData')
    end
%end

%% start main loop here
switch parallelSwitch
    case 'false'
        if numPerm > 0
            if exist([outDir '/permutationResults'],'file') ~= 2
                mkdir([outDir '/permutationResults'])
            end
            permDiv = divisors(numPerm);
            for i = 1:size(permDiv,2)
                numBins(i) = numPerm/permDiv(i);
            end
            permBinSize = numBins(round(size(numBins,2)/2));
            numPermBins = numPerm/permBinSize;
            matFiles = cell(numPermBins,1);
            for permBin = 1:numPermBins
                outFile = [outDir '/permutationResults/Perms_BinSize' num2str(permBinSize) '_BinNumber' num2str(permBin) '.mat'];
                if exist(outFile,'file') == 2
                    load(outFile)
                else
                    matFiles{permBin,1} = outFile;
                    sigmaPerm = zeros(size(SD,1),permBinSize);
                    for perm = 1:permBinSize
                        disp(['Working on perm number ' num2str(perm) ' from bin ' num2str(permBin)]);
                        waitbar((permBin/numPermBins),h,'Running bins of permutations ...')
                        permData = inData(randperm(size(inData,1)),:);
                        switch projectedData.Switch
                            case 'true'
                                warning('Permutation analysis of projected scores is experimental and most likely not accurate ... ')
                                permProjectedScore = projectPCA(origData,'',permData,coeff,'false');
                                permProjectedS = (origData/(permProjectedScore*coeff(:,1:permNumComp)'));
                                permProjectedSD = diag(permProjectedS);
                                permS = permProjectedSD;
                            case 'false'
                                switch centerSwitch
                                    case 'true'
                                        switch pcaType
                                            case 'normal'
                                                permDataCentered = bsxfun(@minus,permData,mean(inData));
                                                [permU, permS, permV] = svd(permDataCentered,'econ');
                                            case 'probabilistic'
                                                permDataCentered = bsxfun(@minus,permData,mean(inData,'omitnan'));
                                                try
                                                    [permU, permS, permV, ~, ~, ~, ~, ~, ~] = ppca_with_svd(permDataCentered,probAttempt);
                                                catch
                                                    disp('Error with PPCA ... data may be rank defficient and inaccurate ... retrying ')
                                                end
                                        end

                                    case 'false'
                                        switch pcaType
                                            case 'normal'
                                                [permU, permS, permV] = svd(permData,'econ');
                                            case 'probabilistic'
                                                [permU, permS, permV, ~, ~, ~, ~, ~, ~] = ppca_with_svd(permData,probAttempt);
                                        end
                                end
                        end
                        switch rotateSwitch
                            case 'true'
                                space = V'*permV;
                                [orthoV,~,orthoU]=svd(space);
                                rotateMat = orthoU*orthoV'; %procrustes rotation
                                permVRotated = permV * permS * rotateMat;
                                permSRotated = sqrt(sum(permVRotated.^2));
                                permSRotatedD = diag(diag(permSRotated));
                                sigmaPerm(:,perm) = permSRotatedD;
                            case 'false'
                                disp('You selected not to perform procrustes rotation ...')
                                permSD = diag(permS);
                                sigmaPerm(:,perm) = permSD;
                        end
                    end
                    save(outFile,'sigmaPerm','matFiles','-v7.3')
                end
            end
            waitbar(0,h,'Getting p-values for permutations ...')
            permPval = nonparamSig(matFiles,'sigmaPerm',SD,'distribution');
            save([outDir '/permutations_pvals.mat'],'permPval')
            %close(h)
        end
        newCompOut = [];
        if numBoot > 0
            clear numBins
            if exist([outDir '/bootstrapResults'],'file') ~= 2
                mkdir([outDir '/bootstrapResults'])
            end
            bootDiv = divisors(numBoot);
            waitbar(0,h,'Starting up bootstrap ...');
            for i = 1:size(bootDiv,2)
                numBins(i) = numBoot/bootDiv(i);
            end
            bootBinSize = numBins(round(size(numBins,2)/2));
            numBootBins = numBoot/bootBinSize;
            matFiles = cell(numBootBins,1);
            for bootBin = 1:numBootBins
                outFile = [outDir '/bootstrapResults/Bootstrap_BinSize' num2str(bootBinSize) '_BinNumber' num2str(bootBin) '.mat'];
                if exist(outFile,'file') == 2
                    load(outFile)
                else
                    matFiles{bootBin,1} = outFile;
                    for boot = 1:bootBinSize
                        disp(['Working on boot number ' num2str(boot) ' from bin ' num2str(bootBin)]);
                        waitbar((bootBin/numBootBins),h,'Running bins of bootstrap ...')
                        switch normalCheckSwitch
                            case 'true'
                                normalCheck = 0;
                                while normalCheck == 0
                                    bootData = bootstrap(inData,1);
                                    normalCheck = kstest(bootData);
                                end
                            case 'false'
                                bootData = bootstrap(inData,1);
                        end
                        %generate the bootstrapped stats
                        switch projectedData.Switch
                            case 'true'
                                switch bootVar
                                    case 'score'
                                        bootScore = projectPCA(origData,'',bootData,coeff,'false');
                                    case 'loadings'
                                        disp('Projection of coefficients based on existing scores not supported ...')
                                        %find a way to get projected
                                        %coefficients
                                    case 'both'
                                        %find a way to get projected
                                        %coefficients
                                        disp('Projection of coefficients based on existing scores not supported running only bootScore...')
                                        bootScore = projectPCA(origData,'',bootData,coeff,'false');
                                end
                            case 'false'
                                switch bootVar
                                    case 'loadings'
                                        switch pcaType
                                            case 'normal'
                                                [bootCoeff,~,~,~,~,~] = pca(bootData,'NumComponents',bootNumComp); %rotations for loadings don't matter since you already have the inputs for rotatefactors
                                            case 'probabilistic'
%                                                 while exist('bootCoeff','var') == 0
                                                    try
                                                        [bootCoeff,~,~,~,~,~] = ppca(bootData,bootNumComp);
                                                    catch ME
                                                        warning(ME.message)
                                                        warning('PPCA will be run with fewer components than specified')
                                                        newCompIdx = strfind(ME.message,'Please try setting k to an integer smaller than ');
                                                        newComp = str2double(ME.message(newCompIdx+48:end-1)) - 1;
                                                        newCompOut = [newCompOut; newComp boot bootBin];
%                                                         if newComp < 6
%                                                             warning('Not enough components')
%                                                         end
                                                        [bootCoeff,~,~,~,~,~] = ppca(bootData,newComp);
                                                    end
%                                                 end
                                        end
                                    case 'scores'
                                        switch pcaType
                                            case 'normal'
                                                [tmpCoeff,bootScore,~,~,~,~] = pca(bootData,'NumComponents',bootNumComp);
                                            case 'probabilistic'
                                                [tmpCoeff,bootScore,~,~,~,~] = ppca(bootData,bootNumComp);
                                        end
                                    case 'both'
                                        switch pcaType
                                            case 'normal'
                                                [bootCoeff,bootScore,~,~,~,~] = pca(bootData,'NumComponents',bootNumComp); %rotations for loadings don't matter since you already have the inputs for rotatefactors
                                            case 'probabilistic'
                                                [bootCoeff,bootScore,~,~,~,~] = ppca(bootData,bootNumComp);
                                        end
                                end
                        end
                        %rotate data if necessary
                        switch rotateSwitch
                            case 'true'
                                switch bootVar
                                    case 'loadings'
                                        [bootCoeff,~] = rotatefactors(bootCoeff,'Method','procrustes','Type','orthogonal','Target',coeff(:,1:size(bootCoeff,2))); %THESE ARE LOADINGS
                                    case 'scores'
                                        [~,rotationMat] = rotatefactors(tmpCoeff,'Method','procrustes','Type','orthogonal','Target',coeff(:,1:size(tmpCoeff,2)));
                                        bootScore = bootScore*rotationMat;
                                    case 'both'
                                        [bootCoeff,rotationMat] = rotatefactors(bootCoeff,'Method','procrustes','Type','orthogonal','Target',coeff(:,1:size(bootCoeff,2))); %THESE ARE LOADINGS
                                        bootScore = bootScore*rotationMat; % THESE ARE SCORES
                                end
                            case 'false'
                                warning('You have chosen not to rotate your scores (by aliging bootstrapped loadings and real loadings with orthogonal procrustes.)')
                        end
                        %okay now write it out
                        switch bootVar
                            case 'loadings'
                                if boot > 1
                                    dataBoot = horzcat(dataBoot,bootCoeff);
                                else
                                    dataBoot = bootCoeff;
                                end
                            case 'scores'
                                if boot > 1
                                    dataBoot = horzcat(dataBoot,bootScore);
                                else
                                    dataBoot = bootScore;
                                end
                            case 'both'
                                if boot > 1
                                    dataBoot = horzcat(dataBoot,bootCoeff);
                                    dataBoot2 = horzcat(dataBoot2,bootScore);
                                else
                                    dataBoot = bootCoeff;
                                    dataBoot2 = bootScore;
                                end
                        end
                    end
                    switch bootVar
                        case 'loadings'
                           % if exist('dataBoot','var') ~= 0
                           if exist('newCompOut','var') ~= 0
                                save(outFile,'dataBoot','matFiles','newCompOut','-v7.3')
                           else
                               save(outFile,'dataBoot','matFiles','-v7.3')
                           end
                           % else
                           %     disp('No data to save...moving onto the next bin')
                           % end
                        case 'scores'
                            save(outFile,'dataBoot','matFiles','-v7.3')
                        case 'both'
                            save(outFile,'dataBoot','dataBoot2','matFiles','-v7.3')
                    end
                end
            end
        end
%         try
%             waitbar(0,h,'Getting p-values for bootstrap ...')
%         catch
%             disp('Getting p-values for bootstrap ... waitbar is broken since you loaded previous data')
%         end
        switch projectedData.Switch
            case 'true'
                switch bootVar
                    case 'loadings'
                        disp('Projection of coefficients based on existing scores not supported ...')
                    case 'scores'
                        waitbar(0,h,'Getting p-values for bootstrap ...')
                        bootPval = nonparamSig(matFiles,'dataBoot',projectedScore,'bootstrap error');
                        save([outDir '/bootstrap_projected_pvals.mat'],'bootPval')
                    case 'both'
                        disp('Projection of coefficients based on existing scores not supported ...')
                        bootPval = nonparamSig(matFiles,'dataBoot',projectedScore,'bootstrap error');
                        save([outDir '/bootstrap_projected_pvals.mat'],'bootPval')
                end
            case 'false'
                switch bootVar
                    case 'loadings'
                        waitbar(0,h,'Getting p-values for bootstrap ...')
                        bootPval = nonparamSig({matFiles{1}},'dataBoot',coeff,'bootstrap error');
                    case 'scores'
                        waitbar(0,h,'Getting p-values for bootstrap ...')
                        bootPval = nonparamSig(matFiles,'dataBoot',score,'bootstrap error');
                    case 'both'
                        waitbar(0,h,'Getting p-values for bootstrap ...')
                        bootPval = nonparamSig(matFiles,'dataBoot',coeff,'bootstrap error');
                        save([outDir '/bootstrap_coeff_pvals.mat'],'bootPval')
                        bootPval2 = nonparamSig(matFiles,'dataBoot2',score,'bootstrap error');
                        save([outDir '/bootstrap_scores_pvals.mat'],'bootPval2')
                end
        end
    case 'true'
        warning('Currently parallel function is not supported due to the way variables are written out')
        if numPerm > 0
            if exist([outDir '/permutationResults'],'file') ~= 2
                mkdir([outDir '/permutationResults'])
            end
            permDiv = divisors(numPerm);
            for i = 1:size(permDiv,2)
                numBins(i) = numPerm/permDiv(i);
            end
            permBinSize = numBins(round(size(numBins,2)/2));
            numPermBins = numPerm/permBinSize;
            matFiles = cell(numPermBins,1);
            for permBin = 1:numPermBins
                outFile = [outDir '/permutationResults/Perms_BinSize' num2str(permBinSize) '_BinNumber' num2str(permBin) '.mat'];
                if exist(outFile,'file') == 2
                    load(outFile)
                else
                    matFiles{permBin,1} = outFile;
                    sigmaPerm = zeros(size(SD,1),permBinSize);
                    parfor perm = 1:permBinSize
                        disp(['Working on perm number ' num2str(perm) ' from bin ' num2str(permBin)]);
                        waitbar((permBin/numPermBins),h,'Running bins of permutations ...')
                        permData = inData(randperm(size(inData,1)),:);
                        switch projectedData.Switch
                            case 'true'
                                warning('Permutation analysis of projected scores is not experimental and most likely not waccurate ... ')
                                permProjectedScore = projectPCA(origData,'',permData,coeff,'false');
                                permProjectedS = (origData/(permProjectedScore*coeff(:,1:permNumComp)'));
                                permProjectedSD = diag(permProjectedS);
                                permS = permProjectedSD;
                            case 'false'
                                switch centerSwitch
                                    case 'true'
                                        permDataCentered = bsxfun(@minus,permData,mean(inData));
                                        [permU, permS, permV] = svd(permDataCentered,'econ');
                                    case 'false'
                                        [permU, permS, permV] = svd(permData,'econ');
                                end
                        end
                        switch rotateSwitch
                            case 'true'
                                space = V'*permV;
                                [orthoV,~,orthoU]=svd(space);
                                rotateMat = orthoU*orthoV'; %procrustes rotation
                                permVRotated = permV * permS * rotateMat;
                                permSRotated = sqrt(sum(permVRotated.^2));
                                permSRotatedD = diag(diag(permSRotated));
                                sigmaPerm(:,perm) = permSRotatedD;
                            case 'false'
                                disp('You selected not to perform procrustes rotation ...')
                                permSD = diag(permS);
                                sigmaPerm(:,perm) = permSD;
                        end
                    end
                    save(outFile,'sigmaPerm','matFiles','-v7.3')
                end
            end
            waitbar(0,h,'Getting p-values for permutations ...')
            permPval = nonparamSig(matFiles,'sigmaPerm',SD,'no');
            save([outDir '/permutations_pvals.mat'],permPval)
            close(h)
        end
        if numBoot > 0
            clear numBins
            if exist([outDir '/bootstrapResults'],'file') ~= 2
                mkdir([outDir '/bootstrapResults'])
            end
            bootDiv = divisors(numBoot);
            waitbar(0,h,'Starting up bootstrap ...');
            for i = 1:size(bootDiv,2)
                numBins(i) = numBoot/bootDiv(i);
            end
            bootBinSize = numBins(round(size(numBins,2)/2));
            numBootBins = numBoot/bootBinSize;
            matFiles = cell(numBootBins,1);
            dataBoot = zeros(size(coeff,1),bootNumComp*bootBinSize);
            dataBoot2 = zeros(size(score,1),bootNumComp*bootBinSize);
            
            for bootBin = 1:numBootBins
                outFile = [outDir '/bootstrapResults/Bootstrap_BinSize' num2str(bootBinSize) '_BinNumber' num2str(bootBin) '.mat'];
                if exist(outFile,'file') == 2
                    load(outFile)
                else
                    matFiles{bootBin,1} = outFile;
                    parfor boot = 1:bootBinSize
                        disp(['Working on boot number ' num2str(boot) ' from bin ' num2str(bootBin)]);
                        waitbar((bootBin/numBootBins),h,'Running bins of bootstrap ...')
                        switch normalCheckSwitch
                            case 'true'
                                normalCheck = 0;
                                while normalCheck == 0
                                    bootData = bootstrap(inData,1);
                                    normalCheck = kstest(bootData);
                                end
                            case 'false'
                                bootData = bootstrap(inData,1);
                        end
                        %generate the bootstrapped stats
                        switch projectedData.Switch
                            case 'true'
                                switch bootVar
                                    case 'score'
                                        bootScore = projectPCA(origData,'',bootData,coeff,'false');
                                    case 'loadings'
                                        disp('Projection of coefficients based on existing scores not supported ...')
                                        %find a way to get projected
                                        %coefficients
                                    case 'both'
                                        %find a way to get projected
                                        %coefficients
                                        disp('Projection of coefficients based on existing scores not supported ...')
                                        bootScore = projectPCA(origData,'',bootData,coeff,'false');
                                end
                            case 'false'
                                switch bootVar
                                    case 'loadings'
                                        %[bootCoeff,~,~,~,~,~] = pca(bootData,'NumComponents',bootNumComp); %rotations for loadings don't matter since you already have the inputs for rotatefactors
                                        [bootCoeff,~,~,~,~,~] = ppca(bootData,bootNumComp);
                                    case 'scores'
                                        %[tmpCoeff,bootScore,~,~,~,~] = pca(bootData,'NumComponents',bootNumComp);
                                        [tmpCoeff,bootScore,~,~,~,~] = ppca(bootData,bootNumComp);
                                    case 'both'
                                        %[bootCoeff,bootScore,~,~,~,~] = pca(bootData,'NumComponents',bootNumComp); %rotations for loadings don't matter since you already have the inputs for rotatefactors
                                        [bootCoeff,bootScore,~,~,~,~] = ppca(bootData,bootNumComp);
                                end
                        end
                        %rotate data if necessary
                        switch rotateSwitch
                            case 'true'
                                switch bootVar
                                    case 'loadings'
                                        [bootCoeff,~] = rotatefactors(bootCoeff,'Method','procrustes','Type','orthogonal','Target',coeff); %THESE ARE LOADINGS
                                    case 'scores'
                                        [~,rotationMat] = rotatefactors(tmpCoeff,'Method','procrustes','Type','orthogonal','Target',coeff);
                                        bootScore = bootScore*rotationMat;
                                    case 'both'
                                        [bootCoeff,rotationMat] = rotatefactors(bootCoeff,'Method','procrustes','Type','orthogonal','Target',coeff); %THESE ARE LOADINGS
                                        bootScore = bootScore*rotationMat; % THESE ARE SCORES
                                end
                            case 'false'
                                warning('You have chosen not to rotate your scores (by aliging bootstrapped loadings and real loadings with orthogonal procrustes.)')
                        end
                        %okay now write it out
                        switch bootVar
                            case 'loadings'
                                
                                
                                
                            case 'scores'
                                if boot > 1
                                    dataBoot = horzcat(dataBoot,bootScore);
                                else
                                    dataBoot = bootScore;
                                end
                            case 'both'
                                if boot > 1
                                    dataBoot = horzcat(dataBoot,bootCoeff);
                                    dataBoot2 = horzcat(dataBoot,bootScore);
                                else
                                    dataBoot = bootCoeff;
                                    dataBoot2 = bootScore;
                                end
                        end
                    end
                    switch bootVar
                        case 'loadings'
                            save(outFile,'dataBoot','matFiles','-v7.3')
                        case 'scores'
                            save(outFile,'dataBoot','matFiles','-v7.3')
                        case 'both'
                            save(outFile,'dataBoot','dataBoot2','matFiles','-v7.3')
                    end
                end
            end
        end
        waitbar(0,h,'Getting p-values for bootstrap ...')
        switch projectedData.Switch
            case 'true'
                switch bootVar
                    case 'loadings'
                        disp('Projection of coefficients based on existing scores not supported ...')
                    case 'scores'
                        bootPval = nonparamSig(matFiles,'dataBoot',projectedScore,'bootstrap error');
                        save([outDir '/bootstrap_projected_pvals.mat'],bootPval)
                    case 'both'
                        disp('Projection of coefficients based on existing scores not supported ...')
                        bootPval = nonparamSig(matFiles,'dataBoot',projectedScore,'bootstrap error');
                        save([outDir '/bootstrap_projected_pvals.mat'],bootPval)
                end
            case 'false'
                switch bootVar
                    case 'loadings'
                        bootPval = nonparamSig(matFiles,'dataBoot',coeff,'bootstrap error');
                    case 'scores'
                        bootPval = nonparamSig(matFiles,'dataBoot',score,'bootstrap error');
                    case 'both'
                        bootPval = nonparamSig(matFiles,'dataBoot',coeff,'bootstrap error');
                        save([outDir '/bootstrap_coeff_pvals.mat'],bootPval)
                        bootPval2 = nonparamSig(matFiles,'dataBoot2',score,'bootstrap error');
                        save([outDir '/bootstrap_scores_pvals.mat'],bootPval2)
                end
        end
end
if exist('permPval','var') == 0
    permPval = [];
end
if exist('bootPval','var') == 0
    bootPval = [];
end
if exist('bootPval','var') == 0
    newCompOut = [];
end
close(h)