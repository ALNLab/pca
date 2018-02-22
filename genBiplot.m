function genBiplot(origLoadings,newLoadings,scorePlot)
% This will generate a biplot. Like matlab's biplot function, data will be
% scaled so that it is easier to visualize (code for scaling and chagning
% sign taken from that script). Unlike matlab's function, it will allow you
% to plot two sets of loadings (from two different analyses) in the same
% space. This may not be meaningful, but if you used projectPCA to get a
% new set of scores, you can use those as your inputs to generate a 3d
% scatter plot of original scores and new scores in the same space. The
% script assumes that you are using the space of the first 3 components.
%
%% Required Inputs:
% origLoadings : this can be your original data. It can be scores or
% loadings.
% newLoadings : these are your new data. Can be left empty if you don't
% want to plot two sets of data. 
% scorePlot : THIS ISN'T FINISHED. PLEASE SET THIS TO 'false'. If set to
% 'true' it will also plot scores as 3D lines. 
%
%% Alex Notes:
% - find a way to plot scores, 
% - make scaling an option,
% - make color an option and give some defaults, 
% - also give support for 1 - 3 components
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Alex Teghipco (alex.teghipco@uci.edu) Last update 5/9/17

warning('Currently does not support plotting scores and loadings simultaneously, input can be either')
if isempty(newLoadings) == 1
    %% Create your own biplot!
    % Force each column of the coefficients to have a positive largest element.
    % This tends to put the large var vectors in the top and right halves of
    % the plot.
    [p,d] = size(origLoadings);
    [~,maxind] = max(abs(origLoadings),[],1);
    colsign = sign(coeff(maxind + (0:p:(d-1)*p)));
    coeff_biplot = bsxfun(@times,origLoadings,colsign);
    
    % Scale the scores so they fit on the plot, and change the sign of
    % their coordinates according to the sign convention for the coefs.
    maxCoefLen = sqrt(max(sum(coeff_biplot.^2,2)));
    score_biplot = bsxfun(@times, maxCoefLen.*(origLoadings ./ max(abs(origLoadings(:)))), colsign);
    score_biplot = origLoadings ./ max(abs(origLoadings(:)));
    
    % Plot away
    figure('Color',[1 1 1]);
    scatter3(score_biplot(:,1),score_biplot(:,2),score_biplot(:,3),'MarkerEdgeColor',[0 0.75 0.85],'MarkerFaceColor',[0 0.75 .85])
    hold on
    scoreLimits = [min(min(score_biplot(:,1:3))) max(max(score_biplot(:,1:3)))];
    scoreLimit = max(abs(scoreLimits));
    scoreRange = [-(scoreLimit):scoreLimit];
    plot3(scoreRange,zeros(size(scoreRange)),zeros(size(scoreRange)),'k')
    plot3(zeros(size(scoreRange)),scoreRange,zeros(size(scoreRange)),'k')
    plot3(zeros(size(scoreRange)),zeros(size(scoreRange)),scoreRange,'k')
    xlabel('Principal Component 1')
    ylabel('Principal Component 2')
    zlabel('Principal Component 3')
    switch scorePlot
        case 'true'
            [p,d] = size(origLoadings);
            zeroes = zeros(p,1); nans = NaN(p,1);
            arx = [zeroes origLoadings(:,1) nans]';
            ary = [zeroes origLoadings(:,2) nans]';
            arz = [zeroes origLoadings(:,3) nans]';
            varHndl = [line(arx(1:2,:),ary(1:2,:),arz(1:2,:), 'Color','b', 'LineStyle','-', 'Marker','none'); ...
                line(arx(2:3,:),ary(2:3,:),arz(2:3,:), 'Color','b', 'Marker','.', plotArgs{:}, 'LineStyle','none')];
            
            set(varHndl(1:p),'Tag','varline');
            set(varHndl((p+1):(2*p)),'Tag','varmarker');
            set(varHndl,{'UserData'},num2cell(([1:p 1:p])'));       
        case 'false'
            disp('Not plotting coefficients...')
    end
else
    %% Create your own biplot!
    % Force each column of the coefficients to have a positive largest element.
    % This tends to put the large var vectors in the top and right halves of
    % the plot.
    [p,d] = size(origLoadings);
    [~,maxind] = max(abs(origLoadings),[],1);
    colsign = sign(coeff(maxind + (0:p:(d-1)*p)));
    coeff_biplot = bsxfun(@times,origLoadings,colsign);
    % Scale the scores so they fit on the plot, and change the sign of
    % their coordinates according to the sign convention for the coefs.
    maxCoefLen = sqrt(max(sum(coeff_biplot.^2,2)));
    score_biplot = bsxfun(@times, maxCoefLen.*(origLoadings ./ max(abs(origLoadings(:)))), colsign);
    score_biplot = origLoadings ./ max(abs(origLoadings(:)));
    
    % Scale the second set of scores so they fit on the plot, and change the sign of
    % their coordinates according to the sign convention for the coefs.
    maxCoefLen = sqrt(max(sum(coeff_biplot.^2,2)));
    score_new_biplot = bsxfun(@times, maxCoefLen.*(projectedLoadings ./ max(abs(projectedLoadings(:)))), colsign);
    score_new_biplot = projectedLoadings ./ max(abs(projectedLoadings(:)));
    
    % Plot away
    figure('Color',[1 1 1]);
    scatter3(score_biplot(:,1),score_biplot(:,2),score_biplot(:,3),'MarkerEdgeColor',[0 0.75 0.85],'MarkerFaceColor',[0 0.75 .85])
    hold on
    scatter3(score_new_biplot(:,1),score_new_biplot(:,2),score_new_biplot(:,3),'MarkerEdgeColor',[0.8 0.4 0],'MarkerFaceColor',[0.8 0.4 0])
    allScores = [score_biplot(:,1:3) score_new_biplot(:,1:3)];
    scoreLimits = [min(min(allScores)) max(max(allScores))];
    scoreLimit = max(abs(scoreLimits));
    scoreRange = [-(scoreLimit):scoreLimit];
    plot3(scoreRange,zeros(size(scoreRange)),zeros(size(scoreRange)),'k')
    plot3(zeros(size(scoreRange)),scoreRange,zeros(size(scoreRange)),'k')
    plot3(zeros(size(scoreRange)),zeros(size(scoreRange)),scoreRange,'k')
    xlabel('Principal Component 1')
    ylabel('Principal Component 2')
    zlabel('Principal Component 3')
    
    switch scorePlot
        case 'true'
            [p,d] = size(origLoadings);
            zeroes = zeros(p,1); nans = NaN(p,1);
            arx = [zeroes origLoadings(:,1) nans]';
            ary = [zeroes origLoadings(:,2) nans]';
            arz = [zeroes origLoadings(:,3) nans]';
            varHndl = [line(arx(1:2,:),ary(1:2,:),arz(1:2,:), 'Color','b', 'LineStyle','-', 'Marker','none'); ...
                line(arx(2:3,:),ary(2:3,:),arz(2:3,:), 'Color','b', 'Marker','.', plotArgs{:}, 'LineStyle','none')];
            
            set(varHndl(1:p),'Tag','varline');
            set(varHndl((p+1):(2*p)),'Tag','varmarker');
            set(varHndl,{'UserData'},num2cell(([1:p 1:p])'));
        case 'false'
            disp('Not plotting coefficients...')
    end
end