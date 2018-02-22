function [topLoadingNames,sortedpercentLoad] = loadingsPlot(inData,thresh,rowNames,colorVec,pVals,pThresh,outDir,trim)
% Plot loadings (i.e. coefficients) from PCA. The plot will organize
% loadings in each component by ascending amount of explained component
% variance and plot the % explained variance for each loading. Positive and
% negative loadings are plotted by different colors in colorVec. An
% arbitrary threshold for interpreting components is denoted on the plot as
% a vertical dashed line. A seperate figure will be generated for each
% component. Each figure will be saved. If you have p-values for each
% loading, only significant loadings will be plotted
%
% inData: rows are loadings, columns are seperate components
%
% thresh: default is 0.25
%
% rowNames: cell array of names that correspond to each loading.
%
% colorVec: needs to be two rows. First row is three columns corresponding
% to the positive loading colors. Second row is three columns corresponding
% to negative loading colors. Assumed to be in RGB.
%
% pVals: if empty, all loadings will be plotted. Otherwise, loadings
% below pThresh only.
%
% pThresh: threshold for p-values (e.g. 0.05, etc)


%% Default parameters
if isempty(outDir) == 0 && exist(outDir,'dir') ~= 7
    mkdir(outDir)
end

if isempty(colorVec) == 1
    colorVec = [1 0.4 0.2;0 0.8 0.8]; %default colors are orange for positive loadings and cyan for negative loadings
else
    if colorVec(1,1) > 1
        colorVec = colorVec/255; %if colovec is supplied we expect it to be RGB
    end
end

if isempty(thresh) == 1
    thresh = 0.25; %default threshold is 0.25
end

if isempty(pThresh) == 1
    pThresh = 0;
end

if isempty(trim) == 1
    trim = 'no';
end

%% change this to change scale. If local, scale is selected by component, if global, scale will be selected across components. 
scale = 'local'; 


%% Convert data to percentage and setup p vals if necessary
% Get percentage of variance explained by ea component and setup axis limits
% get the percentage of component's variance explained. Direction of
% loading doesn't matter for explanatory power.
absLoad = abs(inData);
percentLoad = (absLoad./sum(absLoad))*100;
arbitraryThresh = thresh./sum(absLoad)*100; %how much % variance does a loading of 0.25 explain for ea. component? (arbitrary threshold)

% Remove data below p threshold
if isempty(pVals) == 0
    underP = find(pVals > pThresh);
    percentLoad(underP) = nan;
    inData(underP) = nan;
    threshPercent = min(percentLoad); %ADD THIS TO FUNCTIONALITY
end

% Get axis limits across all figures
switch scale
    case 'global'
    plotLimits = [floor(min(min(percentLoad))) ceil(max(max(percentLoad)))];
end

%% start main loop
for comp = 1:size(inData,2)   
    switch trim
        case 'no'
            compData = inData(:,comp);
            compPercentLoad = percentLoad(:,comp);
            compLabels = rowNames; %this is inside for loop because it will need to be manipulated if p-vals are provided
        case 'yes'
            [absVal,absIdx] = sortrows(absLoad(:,comp));
            if absIdx == 1
                [absVal,absIdx] = sortrows(absLoad(:,comp)');
            end
            absLoadComp = absLoad(absIdx(end-99:end),comp);
            compLabels = rowNames(absIdx(end-99:end));
            compData = inData(absIdx(end-99:end),comp);
            compPercentLoad = percentLoad(absIdx(end-99:end),comp);
    end

    % Get axis limits across all figures
    switch scale
        case 'local'
            plotLimits = [floor(min(min(compPercentLoad))) ceil(max(max(compPercentLoad)))];
    end
        
    %to get around scatter legend bug, scatter plot seperately positive and
    %negative loadings
    negLoad = find(compData<0);
    posLoad = find(compData>0);
    
    % sort percentages
    [sortedpercentLoad(:,comp),sortedpercentLoadIndex] = sort(compPercentLoad,'descend');
    topLoadingNames(:,1) = compLabels(sortedpercentLoadIndex);

    %to fix legend bug in scatter, seperate out now pos/neg loadings and
    %get their associated labels
    for i = 1:size(negLoad,1)
        negLoadSorted(i,1) = find(sortedpercentLoadIndex == negLoad(i));
    end
    for i = 1:size(posLoad,1)
        posLoadSorted(i,1) = find(sortedpercentLoadIndex == posLoad(i));
    end
    sortedPosition(:,1) = (1:size(sortedpercentLoadIndex,1));
    %sortedPosition = flipud(sortedPosition); %this will make the loadings
    %descend instead of ascending
    if isempty(posLoad) == 0
        sortedPositionPos = sortedPosition(posLoadSorted);
        sortedPosLabels = compLabels(posLoad);
    end
    if isempty(negLoad) == 0
        sortedPositionNeg = sortedPosition(negLoadSorted);
        sortedNegLabels = compLabels(negLoad);
    end
    
    %change dot size based on some other magnitude value (to be implemented
    %later)
    dotSizePos = 40;
    dotSizeNeg = 40;
    
    % generate the scatter plot and mind whether you have empty
    % positive/negative loadings vectors
    figure('units','pixels','position',[100 100 1160 990],'color','w')
    if isempty(posLoad) == 0
        scatter((sortedpercentLoad(posLoadSorted,comp)),sortedPositionPos,40,colorVec(1,:),'filled')
        text((sortedpercentLoad(posLoadSorted,comp))+0.02,sortedPositionPos+0.05,sortedPosLabels,'FontSize',7)
        leg = legend('Positive Loadings');
        hold on
    end
    if isempty(negLoad) == 0
        scatter((sortedpercentLoad(negLoadSorted,comp)),sortedPositionNeg,40,colorVec(2,:),'filled')
        text((sortedpercentLoad(negLoadSorted,comp))+0.02,sortedPositionNeg+0.05,sortedNegLabels,'FontSize',7)
        if isempty(posLoad) == 0
            leg = legend('Positive Loadings','Negative Loadings');
        else
            leg = legend('Negative Loadings');
        end
    end

    % add an arrow along y-axis pointing down
    axp = get(gca,'Position');
    xs=axp(1);
    ys=axp(2);
    ye=axp(2)+axp(4)+0.05;
    annotation('arrow', [xs xs],[ye-0.018 ys+0.025]);
    
    % format figure 
    set(gca,'YColor','w') %remove y-axis that overlaps with arrow
    ylabel('Loadings Organized by Ascending Amount of Explained Component Variance','Color','k')
    xlabel('Explained Component Variance')
    set(gca,'XLim',plotLimits)
    set(gca,'YLim',[length(sortedpercentLoad(:,comp)) - sum(~isnan(sortedpercentLoad(:,comp))) length(sortedpercentLoad(:,comp))+1.5]) %increase Ylim for 
    %set(gca,'YLim',[0 size(sortedpercentLoad,1)+1.5])
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',10)
    hold on
    
    textY = length(sortedpercentLoad(:,comp)) - sum(~isnan(sortedpercentLoad(:,comp))) + ((length(sortedpercentLoad(:,comp)) - (length(sortedpercentLoad(:,comp)) - sum(~isnan(sortedpercentLoad(:,comp)))))/2); %find midpoint between data that is not nan
    %textY = sortedpercentLoadIndex(textPt);
    textX = arbitraryThresh(:,comp)-0.02*(plotLimits(2)); %shifted slighlty before the line
    l = line([arbitraryThresh(:,comp) arbitraryThresh(:,comp)],[length(sortedpercentLoad(:,comp)) - sum(~isnan(sortedpercentLoad(:,comp))) length(sortedpercentLoad(:,comp))+1.5],'LineStyle','--','Color',[0.4 0.4 0.4 0.4],'LineWidth',5);
    text(textX,textY,['Salience threshold ' num2str(thresh)],'FontSize',10,'Color',[0.4 0.4 0.4 0.4],'Rotation',90)
    
    %l = line([arbitraryThresh(:,comp) arbitraryThresh(:,comp)],[xs size(sortedpercentLoad,1)],'LineStyle','--','Color',[0.4 0.4 0.4 0.4],'LineWidth',5);
    %text(arbitraryThresh(:,comp)-0.02*(plotLimits(2)),size(sortedpercentLoad,1)/2,['Salience threshold ' num2str(thresh)],'FontSize',10,'Color',[0.4 0.4 0.4 0.4],'Rotation',90)        
    %text(threshPercent+0.2,size(sortedpercentLoad,1)-1,['Salience threshold (' num2str(thresh) ')'],'FontSize',10,'Color',[0.4 0.4 0.4 0.4],'Rotation',90)
    title(['Principal Component ' num2str(comp)],'FontSize',10);   
    leg.String{3} = 'Arbitrary loading signficance threshold';
    
    if exist('threshPercent','var') == 1
       %l2 = line([threshPercent(:,comp) threshPercent(:,comp)],[length(sortedpercentLoad(:,comp)) - sum(~isnan(sortedpercentLoad(:,comp))) length(sortedpercentLoad(:,comp))+1.5],'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',5);
       %text([threshPercent(:,comp)-0.02*(plotLimits(2))],textY,['Salience threshold based on p < ' num2str(pThresh)],'FontSize',10,'Color',[0.2 0.2 0.2],'Rotation',90) 
    end
    
    if isempty(outDir) == 0
        saveas(gcf,[outDir '/Component_' num2str(comp) '_thresh_' num2str(thresh) '.fig'])
        saveas(gcf,[outDir '/Component_' num2str(comp) '_thresh' num2str(thresh) '.eps'],'epsc')
    end
    clear posLoadSorted negLoadSorted sortedPosLabels sortedNegLabels sortedPositionNeg sortedPositionPos negLoad posLoad
    close all
end
