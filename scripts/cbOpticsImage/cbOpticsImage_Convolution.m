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

%% Freeze rng seed and other preamble
rng(1);
markerSizeDownAdjust = 14;

%% Generate some line data
%
% Options:
%   'delta': delta function
%   'shifteddelta': shifted delta function
%   'rand': uniform random noise
%   'chirp: windowed sweep frequency grating
%
% The 'delta' and 'shifteddelta' options are to test that
% everything works as expected for a simple case, while
% the 'chirp' case is to generate an informative figure
% for the book.
data.nLinePixels = 101;
data.smoothFitPoints = 1000;
data.smoothLineX = linspace(0,data.nLinePixels-1,data.smoothFitPoints);
data.lineType = 'chirp';
switch (data.lineType)
    case 'delta'
        data.lineData = zeros(1,data.nLinePixels);
        data.lineData(round(data.nLinePixels/2)) = 1;
    case 'shifteddelta'
        data.lineData = zeros(1,data.nLinePixels);
        data.lineData(round(data.nLinePixels/2)+5) = 1;
    case 'rand'   
        data.lineData = rand(1,data.nLinePixels);
    case 'chirp';
        data.chirpLowFreq = 0;
        data.chirpHighFreq = 0.4;
        data.windowMean = round(data.nLinePixels/2);
        data.windowSd = round(data.nLinePixels/5);
        chirpData = 0.5*chirp(0:data.nLinePixels-1,data.chirpLowFreq,data.nLinePixels,data.chirpHighFreq) + 0.5;
        windowData = normpdf(0:data.nLinePixels-1,data.windowMean,data.windowSd);
        windowData = windowData/max(windowData(:));
        data.lineData = chirpData.*windowData;
end

% Smooth fit to line for plotting
lineFit = fit((0:data.nLinePixels-1)',data.lineData','pchipinterp' );
data.lineSmooth = lineFit(data.smoothLineX')';

%% Generate a line psf.
% 
% Options:
%   'delta': delta function
%   'shifteddelta': shifted delta function
%   'gamma': uniform random noise
%
% The 'delta' and 'shifteddelta' options are to test that
% everything works as expected for a simple case, while
% the 'gamma' case is to generate an informative figure
% for the book.
data.psfType = 'gamma';
data.nLinePsfPixels = 17;
data.linePsfLowPixels = -8;
data.linePsfHighPixels = data.linePsfLowPixels+data.nLinePsfPixels-1;
data.smoothPsfFitPoints = 1000;
data.smoothPsfX = linspace(data.linePsfLowPixels,data.linePsfHighPixels,data.smoothPsfFitPoints);
switch (data.psfType)
    case 'delta'
        data.linePsf = zeros(1,data.nLinePsfPixels);
        data.linePsf(round(data.nLinePsfPixels/2)) = 1;
    case 'shifteddelta'
        data.linePsf = zeros(1,data.nLinePsfPixels);
        data.linePsf(round(data.nLinePsfPixels/2)-5) = 1;
    case 'gamma'
        % Pdf of a gamma function.
        % The various parameters allow one to shift and scale along the x-axis, as
        % well as change the underlying parameters of the gamma itself.
        %
        % The particular parameter choices were made by eye to produce a psf that
        % would work well for the example.  The resulting function is scaled to
        % have a unity area.
        data.linePsfGammaA = 2;
        data.linePsfGammaB = 2.6;
        data.linePsfShrink = 0.5;
        data.linePsfGammaShift = -8;
        data.linePsf = gampdf((data.linePsfLowPixels:data.linePsfHighPixels)/data.linePsfShrink-data.linePsfGammaShift, ...
            data.linePsfGammaA,data.linePsfGammaB);
        data.linePsf = data.linePsf/sum(data.linePsf(:));
end
linePsfFit = fit((data.linePsfLowPixels:data.linePsfHighPixels)',data.linePsf','pchipinterp' );
data.linePsfSmooth = linePsfFit(data.smoothPsfX')';

%% Convolve the signal with the psf
data.blurredLineData = conv2(data.lineData,data.linePsf,'same');
blurredLineFit = fit((0:data.nLinePixels-1)',data.blurredLineData','pchipinterp' );
data.blurredLineSmooth = blurredLineFit(data.smoothLineX')';

%% Now do it in the frequency domain
%
% Take fft of the signal
data.lineFFT = fftshift(fft(data.lineData));

% Make a padded version of the PSF and take the fft of that
psfPadded = zeros(size(data.lineData));
smallCoords = data.linePsfLowPixels:data.linePsfHighPixels;
fullCoords = -round(data.nLinePixels/2):round(data.nLinePixels/2);
for i = 1:length(smallCoords)
    index = find(smallCoords(i) == fullCoords);
    if (~isempty(index))
        psfPadded(index) = data.linePsf(i);
    end
end
data.psfFFT = fftshift(fft(psfPadded));

% Take produce and compute ifft
data.blurredLineFFT = data.lineFFT .* data.psfFFT;
data.fftBlurredLine = ifft(fftshift(data.blurredLineFFT));

%% Make a figure of the line signal, psf, and blurred version.
if (runTimeParams.generatePlots)
    % The line signal
    [lineFigA,figParams] = cbFigInit([],'thinrect');
    figParams.xLimLow = -1;
    figParams.xLimHigh = data.nLinePixels+1;
    figParams.xTicks = [0 25 50 75 100];
    figParams.xTickLabels = {'^{ }0_{ }' '^{ }25_{ }' '^{ }50_{ }' '^{ }75_{ }' '^{ }100_{ }'};
    figParams.yLimLow = 0;
    figParams.yLimHigh = 1;
    figParams.yTicks = [0 0.5 1];
    figParams.yTickLabels = {' 0.0 ' ' 0.5 ' ' 1.0 '};
    
    plot(data.smoothLineX,data.lineSmooth,'r','LineWidth',figParams.lineWidth);
    plot(0:data.nLinePixels-1,data.lineData,'ro','MarkerFaceColor','r','MarkerSize',figParams.markerSize-markerSizeDownAdjust);
    xlabel('Position','FontSize',figParams.labelFontSize);
    ylabel('Image Irrradiance (arbitary units)','FontSize',figParams.labelFontSize);
    title('Input Signal','FontSize',figParams.titleFontSize);
    cbFigAxisSet(lineFigA,figParams);
    
    FigureSave(fullfile(outputDir,[mfilename '_LinePlotA']),lineFigA,figParams.figType);
    
    % The line psf
    [lineFigB,figParams] = cbFigInit([],'thinrect');
    figParams.xLimLow = -8;
    figParams.xLimHigh = 8;
    figParams.xTicks = [-16 -12 -8 -4 0 4 8 12 16];
    figParams.xTickLabels = {'^{ }-8_{ }'  '^{ }-4_{ }' '^{ }0_{ }' '^{ }4_{ }' '^{ }8_{ }'};
    figParams.yLimLow = 0;
    figParams.yLimHigh = 0.5;
    figParams.yTicks = [0 0.25 0.5];
    figParams.yTickLabels = {' 0.00 ' ' 0.25 ' ' 0.50 '};
    
    plot(data.smoothPsfX,data.linePsfSmooth,'r','LineWidth',figParams.lineWidth);
    plot(data.linePsfLowPixels:data.linePsfHighPixels,data.linePsf,'ro','MarkerFaceColor','r','MarkerSize',figParams.markerSize-markerSizeDownAdjust);
    xlabel('Position','FontSize',figParams.labelFontSize);
    ylabel('Point Spread Function','FontSize',figParams.labelFontSize);
    title('Point Spread Function','FontSize',figParams.titleFontSize);
    cbFigAxisSet(lineFigB,figParams);
    
    FigureSave(fullfile(outputDir,[mfilename '_LinePlotB']),lineFigA,figParams.figType);
    
    % The blurred line signal
    [lineFigC,figParams] = cbFigInit([],'thinrect');
    figParams.xLimLow = -1;
    figParams.xLimHigh = data.nLinePixels+1;
    figParams.xTicks = [0 25 50 75 100];
    figParams.xTickLabels = {'^{ }0_{ }' '^{ }25_{ }' '^{ }50_{ }' '^{ }75_{ }' '^{ }100_{ }'};
    figParams.yLimLow = 0;
    figParams.yLimHigh = 1;
    figParams.yTicks = [0 0.5 1];
    figParams.yTickLabels = {' 0.0 ' ' 0.5 ' ' 1.0 '};
    
    plot(data.smoothLineX,data.blurredLineSmooth,'r','LineWidth',figParams.lineWidth);
    plot(0:data.nLinePixels-1,data.blurredLineData,'ro','MarkerFaceColor','r','MarkerSize',figParams.markerSize-markerSizeDownAdjust);
    xlabel('Position','FontSize',figParams.labelFontSize);
    ylabel('Image Irrradiance (arbitary units)','FontSize',figParams.labelFontSize);
    title('Blurred Signal','FontSize',figParams.titleFontSize);
    cbFigAxisSet(lineFigC,figParams);
    
    FigureSave(fullfile(outputDir,[mfilename '_LinePlotC']),lineFigC,figParams.figType);
end



%% Save validation data
%UnitTest.validationData('validateDataStruct', data);

end


