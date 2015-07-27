function varargout = cbOpticsImage_Diffraction(varargin)
%
% Compute and show slices through the diffraction Airy disk.
%
% (c) David Brainard and Andrew Stockman, 2015

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Hello
UnitTest.validationRecord('SIMPLE_MESSAGE', sprintf('%s',mfilename));

%% Set parameters
calcParams.eyeDiameterMm = 24;
calcParams.wavelengthNm = 550;
calcParams.eqCriterionPSFFraction = 0.8;
calcParams.varyPupilDiametersMm = [1 3];
nPupilDiameters = length(calcParams.varyPupilDiametersMm);

%% Compute diffraction blur for each pupil size
% The Psychtoolbox routine AiryPattern does the work.

% Set up radii to compute on.  Start by specifying the range in retinal mm
% and then converting to degrees.  A little roundabout, but whatever.
retinalRadiiMm = 0.5;
retinalRadiiDeg = RetinalMMToDegrees(retinalRadiiMm,calcParams.eyeDiameterMm);
retinalRadiiRad = degtorad(retinalRadiiDeg);

% Set up grid matrices, so that we can convert radius to two-dimensional
% image. Although it is probably inefficient to compute on all the radii of
% a square image matrix (as opposed to computing for linear radii and then
% propogating the andswer onto an image), computers are fast enough that we
% don't care.
nPixels = 501;
centerPixel = round(nPixels+1)/2;
radiusMatrixRaw = MakeRadiusMat(nPixels,nPixels,centerPixel,centerPixel)/nPixels;
radiusMatrixDegs = retinalRadiiDeg*radiusMatrixRaw;
radiusLineDegs = [-radiusMatrixDegs(centerPixel,1:centerPixel-1) radiusMatrixDegs(centerPixel,centerPixel:end)];
radiusMatrixRad = retinalRadiiRad*radiusMatrixRaw;

% Do the calculation for each pupil size and normalize volume of PSF to
% unity.  Also extract 1d slice.
for p = 1:length(calcParams.varyPupilDiametersMm)
    pupilDiameterMm = calcParams.varyPupilDiametersMm(p);
    diffractionPSFImage{p} = AiryPattern(radiusMatrixRad,pupilDiameterMm,calcParams.wavelengthNm); 
    diffractionPSFImage{p} = diffractionPSFImage{p}/sum(diffractionPSFImage{p}(:));
    diffractionPSFSlice{p} = diffractionPSFImage{p}(centerPixel,:);
end

%% Plot slices of the diffraction limited psf
% The plot shows a slice through the center of the psf for each pupil size.
%
% The plot works better to compare shapes if we normalize PSFs to max of 1
% rather than to unit volume, but be aware that the height of the volume
% normalized PSF will be different as a function of pupil size.
if (runTimeParams.generatePlots)
    [diffractionSliceFig,diffractionSliceFigParams] = cbFigInit;
    diffractionSliceFigParams.xLimLow = -0.1;
    diffractionSliceFigParams.xLimHigh = 0.1;
    diffractionSliceFigParams.xTicks = [-0.1 0 0.1];
    diffractionSliceFigParams.xTickLabels = {'^{ }-0.1_{ }' '^{ }0_{ }' '^{ }0.1_{ }'};
    diffractionSliceFigParams.yLimLow = 0;
    diffractionSliceFigParams.yLimHigh = 1;
    diffractionSliceFigParams.yTicks = [0.0 0.2 0.4 0.6 0.8 1];
    diffractionSliceFigParams.yTickLabels = {' 0.0 ' ' 0.2 ' ' 0.4 ' ' 0.6 ' ' 0.8 ' ' 1.0 '};
    
    plot(radiusLineDegs,diffractionPSFSlice{1}/max(diffractionPSFSlice{1}),'r','LineWidth',diffractionSliceFigParams.lineWidth);
    plot(radiusLineDegs,diffractionPSFSlice{2}/max(diffractionPSFSlice{2}),'b','LineWidth',diffractionSliceFigParams.lineWidth);

    xlabel('Retinal Radius (mm)','FontSize',diffractionSliceFigParams.labelFontSize);
    ylabel('Point Spread Function','FontSize',diffractionSliceFigParams.labelFontSize);
    title({'Diffraction Air Disk' ; ' '},'FontSize',diffractionSliceFigParams.titleFontSize);
    cbFigAxisSet(diffractionSliceFig,diffractionSliceFigParams);
    legend({sprintf('^{ } Pupil: %0.2f mm ',calcParams.varyPupilDiametersMm(1)) ...
        sprintf('^{ } Pupil: %0.2f mm ',calcParams.varyPupilDiametersMm(2)) ...
        },'Location','NorthEast','FontSize',diffractionSliceFigParams.legendFontSize);
    FigureSave([mfilename '_VaryPupil'],diffractionSliceFig,diffractionSliceFigParams.figType);
end


%% Save validation data
UnitTest.validationData('calcParams', calcParams);
UnitTest.validationData('diffractionPSFSlice', diffractionPSFSlice);

end







