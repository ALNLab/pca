function [pVal] = nonparamSig(matFiles,varNames,realData,componentAnalysis)
% This is a general function for evaluating the significance of some data
% relative to permutations or bootstrapping that you've run. The script
% takes in either one or many .mat files which are assumed to contain your
% permutation matrix or matrices, and compares them to some real data
% contained in realData. You will need to specify the name of the variable
% within the .mat files that contains the permutation matrix in varNames.
% The columns of your permutation matrix are assumed to correspond to
% permutation or bootstrap samples. The script will automatically
% concatenate all of your samples. If the data does not come from component
% analysis, significance is assesed by comparing real data to the
% distribution of samples, or by using the bootstrap error method. 
%
% Alternatively, the script has some component analysis functionality. For
% instance, it can align your bootstrapped samples with realData. Two
% alignment methods may be used: either an automated procrustes rotation
% that requires the statistics toolbox, and assumes your input data are
% loadings, or a manual method detailed in Peres-Neto, Jackson, and Somers
% (2003; published in Ecology) which accounts for axis reflection and
% reordering by correlating randomized/bootstrapped data to the real data,
% and adjusting component identity based on magnitude of relationship, and
% where necessary, flipping the direction of the relationship. This method
% also assumes loadings are being used. The automated method assesses
% significance using the standard distribution method while the manual
% method requires finding loadings that fit a specific pattern. If a
% loading is positive, then significance is determined by looking at the
% propertion of randomized/bootstrapped loadings that are greater than or
% equal to zero and vice-versa for negative laodings. See the
% aforementioned paper for greater detail.

%% Required Inputs:
% matFiles : a cell with each row corresponding to the full filepath to
% some .mat file that contains your permutations or bootstraps. If this is
% left empty, and you have the script uipickfiles.m in your path, the
% uipickfiles gui will be presented for file selection. In the process of
% concatenation, it is assumed that if you are assessing the significance of
% multiple variables (i.e. your realData has several columns), then these
% variables are all represented in the mat files in the same order. In
% other words if realData has 3 variables, it is assumed that the first 3
% columns in each .mat file correspond to those 3 variables, and so on. 
%
% varNames : this is a string referring to the variable name within the
% matFiles that contains your permutation or bootstrap samples. If
% multiple variables from each .mat file need to be concatenated, then this
% is a cell where each row corresponds to one such variable name. Variable
% names must exist in all of your .mat files. 
%
% realData : this is an array of your realData. Columns represent variables
% and rows represent data points. If you are testing significance for
% several variables, then it is assumed that your mat files will contain
% several variables. See matFiles for assumptions about variable
% locations in randomized/bootstrapped data. 
%
% componentAnalysis : this is a switch that gives several options to assess
% significance of data -- if set to 'distribution', the standard method
% will be used of counting the number of values in the
% randomized/bootstrapped data that is greater than the real values. If it
% is set to 'bootstrap error', the real data will be divided by the
% standard deviation of the bootstrapped distribution to generate z-scores
% that are assumed to be two-tail and converted to p-values. If set to
% 'bootstrapped loadings (manual rotation)' the Peres-Neto, Jackson, and
% Somers (2003; published in Ecology) method is used, and if set to
% 'bootstrapped loadings (auto rotation)', the bootstrap error method is
% used AFTER aligning your bootstrapped loadings to the real loadings with
% an ortho procrustes rotation. In the case of manual rotation, signficance
% is assessed by the number of loadings showing the correct trend. If a
% loading is positive, then significance is determined by looking at the
% propertion of randomized/bootstrapped loadings that are greater than or
% equal to zero and vice-versa for negative laodings. See the
% aforementioned paper for greater detail.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Alex Teghipco (alex.teghipco@uci.edu)
% Last update 5/4/17


if isempty(matFiles) == 1
    matFiles = uipickfiles';
end
if ischar(varNames) == 1
    extractVar = cell(size(matFiles,1),1);
    for i = 1:size(extractVar,1)
        extractVar{i,1} = varNames;
    end
end

for i = 1:size(extractVar,1)
    disp([num2str((i/size(extractVar,1))*100) '% of matrix concatenation complete...'])
    vars = extractVar(i,:);
    vars = vars(~cellfun('isempty',vars));
    for j = 1:size(vars,2)
        if j == 1
            inMat = load(matFiles{i});
            tmp = inMat.(vars{j});
            if isfield(inMat,'newCompOut') ~= 0 
                tmp2 = inMat.newCompOut;
            end
        else
            tmp = horzcat(tmp,inMat.(vars{j}));
            if isfield(inMat,'newCompOut') ~= 0
                tmp2 = horzcat(tmp2,inMat.newCompOut);
            end
        end
    end
    if i == 1
        if exist('tmp2','var') ~= 0 
           concatComp = tmp2; 
        end
        concatMat = tmp;
    else
        if exist('tmp2','var') ~= 0
            concatComp = horzcat(concatComp,tmp2);
        end
        concatMat = horzcat(concatMat,tmp);
    end
    clear tmp
end
switch componentAnalysis
    case 'distribution'
        disp('Testing significance...')
        for i = 1:size(realData,1)
            idx = find(realData(i,1) >= concatMat(i,:));
            pVal(i,1) = 1 - (size(idx,2)/size(concatMat,2));
            clear idx
        end
    case 'bootstrapped loadings (manual rotation)'
        disp('Correcting for axis reflection and reordering ...')
        concatMatR = reshape(concatMat,[size(score,1),size(score,2),(size(concatMat,2)/size(score,2))]);
        bootstatCorr = nan(size(concatMatR));
        for i = 1:size(concatMatR,3)
            boot_tmp = concatMatR(:,:,i);
            tmpR = [boot_tmp realData];
            [rMat,~] = corrcoef(tmpR,'rows','pairwise');
            for j = 1:size(realData,2)
                [~,trueIdx(j)] = max(abs(rMat(j,size(boot_tmp,2)+1:size(boot_tmp,2)+3)));
                [trueVal(j),~] = max(rMat(j,size(boot_tmp,2)+1:size(boot_tmp,2)+3));
            end
            invIdx = find(trueVal < 0);
            boot_tmp(:,invIdx) = boot_tmp(:,invIdx)*-1;
            
            [~,row] = sort(trueIdx);
            bootstatCorr(:,:,i) = boot_tmp(:,row);
        end
        disp('Testing significance ...')
        for i = 1:size(realData,1)
            for j = 1:size(realData,2)
                loading = realData(i,j);
                if loading > 0
                    idx = find(bootstatCorr(i,j,:) <= 0);
                else
                    idx = find(bootstatCorr(i,j,:) >= 0);
                end
                pVal(i,j) = 1 - (size(idx,1)/size(bootstatCorr,2));
            end
        end
    case 'bootstrapped loadings (auto rotation)'
        disp('Assuming you have already checked for normality of boostrap distribution ...')
        disp('Running procrustes rotation on loadings  ...')
        concatMatR = reshape(concatMat,[size(realData,1),size(realData,2),(size(concatMat,2)/size(score,2))]);
        bootstatCorr = nan(size(concatMatR));
        for i = 1:size(concatMatR,2)
            concatMatRR(:,:,i) = rotatefactors(concatMatR(:,:,i),'Method','procrustes','Type','orthogonal','Target',realData);
        end
        disp('Testing significance ...')
        bootError = std(concatMatRR,0,3,'omitnan');
        zScore = score./bootError;
        pVal = 2*(1-normcdf(abs(zScore),0,1));
    case 'bootstrap error'
        disp('Assuming you have already checked for normality of boostrap distribution ...')
        disp('Reshaping ... ')
        %tmpLoad
        if exist('concatComp','var') == 0
            concatMatR = reshape(concatMat,[size(realData,1),size(realData,2),(size(concatMat,2)/size(realData,2))]);
            bootError = std(concatMatR,0,3,'omitnan');
            zScore = realData./bootError;
            pVal = 2*(1-normcdf(abs(zScore),0,1));
        else
            comp = concatComp(:,1:3:end);
            compR = reshape(comp,[size(comp,2)*size(comp,1),1]);
            boot = concatComp(:,2:3:end);
            bootR = reshape(boot,[size(boot,2)*size(boot,1),1]);
            bin = concatComp(:,3:3:end);
            binR = reshape(bin,[size(bin,2)*size(bin,1),1]);
            
            for i = 1:size(compR,1)
                compEnd(i,1) = sum(compR(1:i));
                if i == 1
                    compStart(i,1) = 1;
                else
                    compStart(i,1) = compEnd(i-1,1) + 1;
                end
            end
            
            for compi = 1:max(compR)
                concatIdx = find(compR >= compi);
                for j = 1:size(concatIdx,1)
                    compiCols(j,1) = compStart(concatIdx(j))+compi-1;
                end
                bootError(:,compi) = std(concatMat(:,compiCols),0,2,'omitnan');
                clear compiCols               
            end
            
            zScore = realData(1:size(bootError,2))./bootError;
            pVal = 2*(1-normcdf(abs(zScore),0,1));

        end

end
