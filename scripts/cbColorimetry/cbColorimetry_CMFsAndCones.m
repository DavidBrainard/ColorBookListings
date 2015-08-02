function varargout = cbColorimetry_CMFsAndCones(varargin)
%
% Illustrates connections between color matching functions and cone fundamentals.
%
%
% (c) David Brainard and Andrew Stockman, 2015

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Hello
UnitTest.validationRecord('SIMPLE_MESSAGE', sprintf('%s',mfilename));
outputDir = sprintf('%s_Output',mfilename);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

%% Set parameters

%% Get and plot Stiles-Burch 10-degree color matching functions
load T_stiles10;
wls = SToWls(S_stiles10);
data.wls = wls;
data.S_stiles10 = S_stiles10;
data.T_stiles10 = T_stiles10;
clear S_stiles10 T_stiles10
if (runTimeParams.generatePlots)
    [stilesBurch10Fig,figParams] = cbFigInit;
    figParams.xLimLow = 350;
    figParams.xLimHigh = 750;
    figParams.xTicks = [350 400 450 500 550 600 650 700 750];
    figParams.xTickLabels = {'^{ }350_{ }' '^{ }400_{ }' '^{ }450_{ }' '^{ }500_{ }' ...
        '^{ }550_{ }' '^{ }600_{ }' '^{ }650_{ }' '^{ }700_{ }' '^{ }750_{ }'};
        figParams.yLimLow = -1;
    figParams.yLimHigh = 4;
    figParams.yTicks = [-1 0 1 2 3 4];
    figParams.yTickLabels = {'-1.0 ' ' 0.0 ' ' 1.0 ' ' 2.0 ' ' 3.0 ' ' 4.0 '};
    
    plot(data.wls,data.T_stiles10(1,:)','r','LineWidth',figParams.lineWidth);
    plot(data.wls,data.T_stiles10(2,:)','g','LineWidth',figParams.lineWidth);
    plot(data.wls,data.T_stiles10(3,:)','b','LineWidth',figParams.lineWidth);
    
    xlabel('Wavelength (nm)','FontSize',figParams.labelFontSize);
    ylabel('CMF (energy units)','FontSize',figParams.labelFontSize);
    title('Stiles-Burch 10-degree CMFs','FontSize',figParams.titleFontSize);
    cbFigAxisSet(stilesBurch10Fig,figParams);
 
    % Save the figure
    FigureSave(fullfile(outputDir,[mfilename '_StilesBurch10Cmfs']),stilesBurch10Fig,figParams.figType);
end

%% Load Stockman-Sharpe 10-degree cone fundamentals
load T_cones_ss10
data.T_cones10 = SplineCmf(S_cones_ss10,T_cones_ss10,wls);
data.M_CmfToCones = ((data.T_stiles10')\(data.T_cones10'))';
data.T_cones10_fit = data.M_CmfToCones*data.T_stiles10;
if (runTimeParams.generatePlots)
    [stockmanSharpe10Fig,figParams] = cbFigInit;
    figParams.xLimLow = 350;
    figParams.xLimHigh = 750;
    figParams.xTicks = [350 400 450 500 550 600 650 700 750];
    figParams.xTickLabels = {'^{ }350_{ }' '^{ }400_{ }' '^{ }450_{ }' '^{ }500_{ }' ...
        '^{ }550_{ }' '^{ }600_{ }' '^{ }650_{ }' '^{ }700_{ }' '^{ }750_{ }'};
    figParams.yLimLow = 0;
    figParams.yLimHigh = 1;
    figParams.yTicks = [0 0.5 1];
    figParams.yTickLabels = {' 0.0 ' ' 0.5 ' ' 1.0 '};
    
    plot(data.wls,data.T_cones10(1,:)','r','LineWidth',figParams.lineWidth);
    plot(data.wls,data.T_cones10(2,:)','g','LineWidth',figParams.lineWidth);
    plot(data.wls,data.T_cones10(3,:)','b','LineWidth',figParams.lineWidth);
    
    plot(data.wls,data.T_cones10_fit(1,:)','k:','LineWidth',figParams.lineWidth-1);
    plot(data.wls,data.T_cones10_fit(2,:)','k:','LineWidth',figParams.lineWidth-1);
    plot(data.wls,data.T_cones10_fit(3,:)','k:','LineWidth',figParams.lineWidth-1); 
    
    xlabel('Wavelength (nm)','FontSize',figParams.labelFontSize);
    ylabel('Cone Fundamental (energy units)','FontSize',figParams.labelFontSize);
    title('Stiles-Burch 10-degree CMFs','FontSize',figParams.titleFontSize);
    cbFigAxisSet(stockmanSharpe10Fig,figParams);
    
    % Save the figure
    FigureSave(fullfile(outputDir,[mfilename '_StockmanSharpe10']),stockmanSharpe10Fig,figParams.figType);
end

%% Make a 3D plot of the spectrum locus
S_10nm = [390 10 37];
data.T_stiles10_10nm = SplineCmf(data.wls,data.T_stiles10,S_10nm);
if (runTimeParams.generatePlots)
    [stilesBurch10SpectrumLocusFig,figParams] = cbFigInit;
    figParams.xLimLow = -1;
    figParams.xLimHigh = 4;
    figParams.xTicks = [-1 0 1 2 3 4];
    figParams.xTickLabels = {'-1.0 ' ' 0.0 ' ' 1.0 ' ' 2.0 ' ' 3.0 ' ' 4.0 '};
    figParams.yLimLow = -0.5;
    figParams.yLimHigh = 1.5;
    figParams.yTicks = [-0.5 0 0.5 1.0 1.5 ];
    figParams.yTickLabels = {'-0.5 ' ' 0.0 ' ' 0.5 ' ' 1.0 ' ' 1.5 '};
    figParams.zLimLow = -0.5;
    figParams.zLimHigh = 1.5;
    figParams.zTicks = [-0.5 0 0.5 1.0 1.5 ];
    figParams.zTickLabels = {'-0.5' ' 0.0 ' ' 0.5 ' ' 1.0 ' ' 1.5 '};
    
    plot3(data.T_stiles10(1,:)',data.T_stiles10(2,:)',data.T_stiles10(3,:)', ...
        'k','LineWidth',figParams.lineWidth);
    plot3(data.T_stiles10_10nm(1,:)',data.T_stiles10_10nm(2,:)',data.T_stiles10_10nm(3,:)', ...
        'ko','MarkerFaceColor','k','MarkerSize',figParams.markerSize-10);
    
    xlabel('r','FontSize',figParams.labelFontSize);
    ylabel('g','FontSize',figParams.labelFontSize);
    zlabel('b','FontSize',figParams.labelFontSize);
    title('Stiles-Burch 10-degree CMFs','FontSize',figParams.titleFontSize);
    cbFigAxisSet(stilesBurch10SpectrumLocusFig,figParams);
    zlim([figParams.zLimLow figParams.zLimHigh]);
    set(gca,'ZTick',figParams.zTicks);
    set(gca,'ZTickLabels',figParams.zTickLabels);
    az = 55; el = 26; view(az,el);
    grid on
 
    % Save the figure
    FigureSave(fullfile(outputDir,[mfilename '_StilesBurch10SpectrumLocus']),stilesBurch10SpectrumLocusFig,figParams.figType);
end

%% Save validation data
UnitTest.validationData('validateDataStruct', data);

end

% Legend, with tweak to make lines long enough so that dash shows.
% Note the extra spaces that preface the actual legend text. Ugh.
% [~,legendChildObjs] = legend({['^{ }' figParams.legendExtraSpaceStr '  Linear '],[ '^{ }' figParams.legendExtraSpaceStr '  Model Eye Based ']},...
%     'Location','NorthWest','FontSize',figParams.legendFontSize);
% lineObjs = findobj(legendChildObjs, 'Type', 'line');
% xCoords = get(lineObjs, 'XData') ;
% for lineIdx = 1:length(xCoords)
%     if (length(xCoords{lineIdx}) ~= 2), continue; end
%     set(lineObjs(lineIdx), 'XData', xCoords{lineIdx} + [0 figParams.legendLineTweak])
% end


