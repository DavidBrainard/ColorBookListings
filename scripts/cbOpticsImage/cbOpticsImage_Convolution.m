function varargout = cbOpticsImage_Convolution(varargin)
%
% Produce figures to illustrate convolution.
%
% Does both a one dimensional and two dimensional example, each single
% channel.
%
% (c) David Brainard and Andrew Stockman, 2016

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
close all; drawnow;

%% Load in an image
%
% Grab a hyperspectral scene from our archiva server
% and go from there.

% Load a scene and extract a normalized grayscale image
client = RdtClient('isetbio');
remotePath = 'resources/scenes/hyperspectral/manchester_database';
client.crp(remotePath);
[theScene, theSceneArtifact] = client.readArtifact('scene7', 'type', 'mat');
imageData = theScene.scene.data.photons(:,:,15);

% Quick and dirty tone map
upFactor = 1.5;
imageData = imageData/(upFactor*max(imageData(:)));
imageData = sqrt(imageData);
imageDdata = imageData.^2;
imageData = imageData/max(imageData(:));

% Extract line data
nPixels = size(imageData,1);
lineData = imageData(round(nPixels/2),:);

% Take a look at the image, sanity check
imshow(imageData);
figure

%% Make a figure of the spectrum and its fit for the box
if (runTimeParams.generatePlots)
    [spectralFig,figParams] = cbFigInit;
    figParams.xLimLow = 1;
    figParams.xLimHigh = nPixels;
    figParams.xTicks = [350 400 450 500 550 600 650 700 750 800];
    figParams.xTickLabels = {'^{ }350_{ }' '^{ }400_{ }' '^{ }450_{ }' '^{ }500_{ }' ...
        '^{ }550_{ }' '^{ }600_{ }' '^{ }650_{ }' '^{ }700_{ }' '^{ }750_{ }' '^{ }800_{ }'};
    figParams.yLimLow = 0;
    figParams.yLimHigh = 1.5;
    figParams.yTicks = [0 0.5 1 1.5];
    figParams.yTickLabels = {' 0.0 ' ' 0.5 ' ' 1.0 ' ' 1.5 '};
    
    plot(data.theWls,data.theSpdIrradiance,'r','LineWidth',figParams.lineWidth);
    plot(data.theWls,data.theSpdSynthesized,'k:','LineWidth',figParams.lineWidth-1);
    
    xlabel('Wavelength (nm)','FontSize',figParams.labelFontSize);
    ylabel('Irradiance (Watts/[m2-nm])','FontSize',figParams.labelFontSize);
    cbFigAxisSet(spectralFig,figParams);
    [~,legendChildObjs] = legend({['^{ }' figParams.legendExtraSpaceStr '  Spectrum  '],[ '^{ }' figParams.legendExtraSpaceStr '  Synthesized Spectrum']}, ...
        'Location','NorthEast','FontSize',figParams.legendFontSize);
    lineObjs = findobj(legendChildObjs, 'Type', 'line');
    xCoords = get(lineObjs, 'XData') ;
    for lineIdx = 1:length(xCoords)
        if (length(xCoords{lineIdx}) ~= 2), continue; end
        set(lineObjs(lineIdx), 'XData', xCoords{lineIdx} + [0 figParams.legendLineTweak])
    end
    FigureSave(fullfile(outputDir,[mfilename '_Spectrum']),spectralFig,figParams.figType);
end

%% Figure of the scaled narrowband lights
if (runTimeParams.generatePlots)
    [narrowbandFig,figParams] = cbFigInit;
    figParams.xLimLow = 350;
    figParams.xLimHigh = 800;
    figParams.xTicks = [350 400 450 500 550 600 650 700 750 800];
    figParams.xTickLabels = {'^{ }350_{ }' '^{ }400_{ }' '^{ }450_{ }' '^{ }500_{ }' ...
        '^{ }550_{ }' '^{ }600_{ }' '^{ }650_{ }' '^{ }700_{ }' '^{ }750_{ }' '^{ }800_{ }'};
    figParams.yLimLow = 0;
    figParams.yLimHigh = 1.5;
    figParams.yTicks = [0 0.5 1 1.5];
    figParams.yTickLabels = {' 0.0 ' ' 0.5 ' ' 1.0 ' ' 1.5 '};
    
    plot(data.theWls,data.scaledB,'LineWidth',figParams.lineWidth);
    
    xlabel('Wavelength (nm)','FontSize',figParams.labelFontSize);
    ylabel('Irradiance (Watts/[m2-nm])','FontSize',figParams.labelFontSize);
    cbFigAxisSet(spectralFig,figParams);
    %legend({'Linear', 'Model Eye Based'},'Location','NorthWest','FontSize',figParams.legendFontSize);
    FigureSave(fullfile(outputDir,[mfilename '_NarrowbandSpectra']),spectralFig,figParams.figType);
end

%% Save validation data
UnitTest.validationData('validateDataStruct', data);

end


