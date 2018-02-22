function [projected] = projectPCA(origData,origScores,newData,origCoeffs,plotSwitch)
% This function will generate scores for an entirely new dataset based on
% loadings from an original dataset, thereby putting the new dataset into
% the pca space from the first dataset. It can also generate a biplot of
% scores taken from the original dataset, and scores taken from the new
% dataset. 
%
%% Required Inputs:
% origData : this is your original data. Columns are taken to be variables,
% rows as observations.
% origScores : these are your original scores, with columns representing
% components. This can be left empty if you are not planning on generating
% a biplot.
% newData : this is your new data. Columns are taken to be variables, and
% rows to be observations.
% origCoeffs : these are the loadings from your original data. Columns are
% taken to be components.
% plotSwitch : if set to 'true, a biplot will be generated overlaying the
% original scores and the new scores. Otherwise, set to 'false'.
%
%% Alex Notes
% - MAKE GENPLOT PLOT COEFFICIENTS AS WELL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Alex Teghipco (alex.teghipco@uci.edu)
% Last update 5/4/17


centered = bsxfun(@minus, newData, mean(origData));
projected = centered*origCoeffs; %original method, this gives scores

switch plotSwitch
    case 'true'
        warning('Plotting coefficients not yet supported.')
        genBiplot(origCoeffs,origScores,newLoadings,'true')
    case 'false'
        disp('Finished projecting new data matrix into old pca space ...')
end
