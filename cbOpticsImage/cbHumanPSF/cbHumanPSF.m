function cbHumanPSF
% cbPinholeOpticsBlur
%
% Compute PSF for human eye.
%
% Requires: Psychophysics Toolbox, isetbio
%
% (c) David Brainard and Andrew Stockman, 2014

%% Clear
clear; close all;

%% Set parameters
eyeLengthMm = 17;
wavelengthNm = 550;
distanceToSourceMm = 2000;
eqCriterionPSFFraction = 0.8;
pupilDiametersMm = [0.05 0.1 0.2];
nPupilDiameters = length(pupilDiametersMm);

%% Compute diffraction blur for each pupil size
% The Psychtoolbox routine AiryPattern does the work.

% Set up radii to compute on.  Start by specifying the range in retinal mm
% and then converting to degrees.
retinalRadiiMm = 0.5;
retinalRadiiDeg = RetinalMMToDegrees(retinalRadiiMm,eyeLengthMm);
retinalRadiiRad = degtorad(retinalRadiiDeg);

% Set up grid matrices, so that we can convert radius to two-dimensional
% image. Although it is probably inefficient to compute on all the radii of
% a square image matrix (as opposed to computing for linear radii and then
% propogating the andser onto an image), computers are fast enough that we
% don't care.
nPixels = 501;
centerPixel = round(nPixels+1)/2;
radiusMatrixRaw = MakeRadiusMat(nPixels,nPixels,centerPixel,centerPixel)/nPixels;
radiusMatrixDegs = retinalRadiiDeg*radiusMatrixRaw;
radiusMatrixRad = retinalRadiiRad*radiusMatrixRaw;
radiusMatrixMm = retinalRadiiMm*radiusMatrixRaw;
radiusLineMm = radiusMatrixMm(centerPixel,centerPixel:end);

% Do the calculation for each pupil size and normalize volume of PSF to
% unity.  Also extract 1d slice.
for p = 1:length(pupilDiametersMm)
    pupilDiameterMm = pupilDiametersMm(p);
    diffractionPSFImage{p} = AiryPattern(radiusMatrixRad,pupilDiameterMm,wavelengthNm); 
    diffractionPSFImage{p} = diffractionPSFImage{p}/sum(diffractionPSFImage{p}(:));
    diffractionPSFSlice{p} = diffractionPSFImage{p}(centerPixel,centerPixel:end);
end

%% Compute equivalent blur circle
% For comparison with geometric blur, it is convenient to characterize the
% diffraction limited PSF by an equivalent blur circle.  We do this by
% finding the radius that contains a criterion fraction of the pupil
% volume, and calling that the equivlent cirular psf. This is a rough and
% ready approximation, but we find it conceptually convenient as a summary
% of the size of the PSF.
radiiMm = unique(radiusMatrixMm(:));
for p = 1:length(pupilDiametersMm)
    for i = 2:length(radiiMm)
        index = find(radiusMatrixMm <= radiiMm(i));
        volume(i) = sum(diffractionPSFImage{p}(index));
        if (volume(i) > eqCriterionPSFFraction)
            lambda = (eqCriterionPSFFraction-volume(i-1))/(volume(i)-volume(i-1));
            eqDiffractionBlurCircleDiameterMm(p) = (1-lambda)*radiiMm(i-1) + lambda*radiiMm(i);
            eqDiffractionBlurCircleDiameterDegs(p) = RetinalMMToDegrees(eqDiffractionBlurCircleDiameterMm(p),eyeLengthMm);
            break;
        end
    end
    
    % Compute circular psfs at the equivalent diameters
    %
    % Build the image
    eqDiffractionPSFImageMm{p} = ones(size(radiusMatrixMm));
    index = find(radiusMatrixMm > eqDiffractionBlurCircleDiameterMm(p));
    eqDiffractionPSFImageMm{p}(index) = 0;
    
    % Normalize volume and extract slice
    eqDiffractionPSFImageMm{p} = eqDiffractionPSFImageMm{p}/sum(eqDiffractionPSFImageMm{p}(:));
    eqDiffractionPSFSlice{p} = eqDiffractionPSFImageMm{p}(centerPixel,centerPixel:end);
    
    % Print summary of this calculation
    fprintf('Pupil size %0.2f mm, diffraction equiv blur cicle (%d%% volume) %0.3f mm, %0.3f deg\n',...
        pupilDiametersMm(p),round(100*eqCriterionPSFFraction),eqDiffractionBlurCircleDiameterMm(p),eqDiffractionBlurCircleDiameterDegs(p));
end
fprintf('\n');

%% Plot a slice of the diffraction limited psf
% The plot shows a slice through the center of the psf for two pupil sizes
% (the smallest and largest that we compute for.).
%
% The plot works better to compare shapes if we normalize PSFs to max of 1
% rather than to unit volume, but be aware that the height of the volume
% normalized PSF will be different as a function of pupil size.
%
% The plot also shows radius of equivalent blur circle as dashed vertical
% lines.
[diffractionSliceFig,diffractionSliceFigParams] = cbFigInit;
diffractionSliceFigParams.xLimLow = 0;
diffractionSliceFigParams.xLimHigh = 0.4;
diffractionSliceFigParams.xTicks = [0 0.1 0.2 0.3 0.4];
diffractionSliceFigParams.xTickLabels = {};
diffractionSliceFigParams.yLimLow = 0;
diffractionSliceFigParams.yLimHigh = 1;
diffractionSliceFigParams.yTicks = [0.0 0.2 0.4 0.6 0.8 1];
diffractionSliceFigParams.yTickLabels = {};
plot(radiusLineMm,diffractionPSFSlice{1}/max(diffractionPSFSlice{1}),'r','LineWidth',diffractionSliceFigParams.lineWidth);
plot(radiusLineMm,diffractionPSFSlice{end}/max(diffractionPSFSlice{end}),'b','LineWidth',diffractionSliceFigParams.lineWidth);
plot([eqDiffractionBlurCircleDiameterMm(1) eqDiffractionBlurCircleDiameterMm(1)],[0 0.5],'r:','LineWidth',diffractionSliceFigParams.lineWidth-1);
plot([eqDiffractionBlurCircleDiameterMm(end) eqDiffractionBlurCircleDiameterMm(end)],[0 0.5],'b:','LineWidth',diffractionSliceFigParams.lineWidth-1);
xlabel('Retinal Radius (mm)','FontSize',diffractionSliceFigParams.labelFontSize);
ylabel('Point Spread Function','FontSize',diffractionSliceFigParams.labelFontSize);
title('Pinhole Camera - Diffraction Limited Blur','FontSize',diffractionSliceFigParams.titleFontSize);
cbFigAxisSet(diffractionSliceFig,diffractionSliceFigParams);
legend({sprintf('Pupil: %0.2f mm',pupilDiametersMm(1)) sprintf('Pupil: %0.2f mm',pupilDiametersMm(end))},'Location','NorthEast','FontSize',diffractionSliceFigParams.legendFontSize);
FigureSave('PinholeOpticsBlurDiffractionSlice',diffractionSliceFig,diffractionSliceFigParams.figType);

%% Compute geometric blur for a pinhole optics.
% This depends on the distance to the object, and in the limit of a
% infitely distant point source is just the pupil diameter directly.
%
% We think that the distance dependence is also true of diffraction, in the
% sense that using the Airy pattern as the PSF results from some
% approximations that treat the arriving wavefront as planar at the pupil.
%
% In any case, we'll use a distance that is big with respect to the scale
% of the model eye.
for p = 1:length(pupilDiametersMm)
    % Geometric calculation
    geometricBlurCircleDiameterMm(p)  = ((distanceToSourceMm+eyeLengthMm)/distanceToSourceMm)*pupilDiametersMm(p);
    
    % For a really fair comparison with diffraction, should find the
    % equivalent circle diameter, that contains the criterion fraction of
    % the volume.
    eqGeometricBlurCircleDiameterMm(p) = sqrt(eqCriterionPSFFraction)*geometricBlurCircleDiameterMm(p);
    eqGeometricBlurCircleDiameterDegs(p) = RetinalMMToDegrees(eqGeometricBlurCircleDiameterMm(p),eyeLengthMm);
    eqGeometricPSFImageMm{p} = ones(size(radiusMatrixMm));
    index = find(radiusMatrixMm > geometricBlurCircleDiameterMm(p));
    eqGeometricPSFImageMm{p}(index) = 0;
    
    % Normalize volume and extract slice
    eqGeometricPSFImageMm{p} = eqGeometricPSFImageMm{p}/sum(eqGeometricPSFImageMm{p}(:));
    eqGeometricPSFSlice{p} = eqGeometricPSFImageMm{p}(centerPixel,centerPixel:end);
    
     % Print summary of this calculation
    fprintf('Pupil size %0.2f mm, geometric equiv blur cicle (%d%% volume) %0.3f mm, %0.3f deg\n',...
        pupilDiametersMm(p),round(100*eqCriterionPSFFraction),eqGeometricBlurCircleDiameterMm(p),eqGeometricBlurCircleDiameterDegs(p));
end
fprintf('\n');

%% Plot a slice through the geometric blur circle
% The plot shows the geometry-limited PSF, which is just a circle. It
% doesn't look quite like a circle because of numerical precision issues in
% the 2D computation of the PSF implemented here, but the basic point is
% clear.
%
% The plot also shows as a dashed line the radius that contains the same
% criterion fraction of the PSF mass as for the diffraction limited
% calculation.  This provides a metric that may be compared to the size of
% the same metric for the diffraction limited PSF.
[geometricSliceFig,geometricSliceFigParams] = cbFigInit;
geometricSliceFigParams.xLimLow = 0;
geometricSliceFigParams.xLimHigh = 0.4;
geometricSliceFigParams.xTicks = [0 0.1 0.2 0.3 0.4];
geometricSliceFigParams.xTickLabels = {};
geometricSliceFigParams.yLimLow = 0;
geometricSliceFigParams.yLimHigh = 1;
geometricSliceFigParams.yTicks = [0.0 0.2 0.4 0.6 0.8 1];
geometricSliceFigParams.yTickLabels = {};
plot(radiusLineMm,eqGeometricPSFSlice{1}/max(eqGeometricPSFSlice{1}),'r','LineWidth',geometricSliceFigParams.lineWidth+1);
plot(radiusLineMm,eqGeometricPSFSlice{end}/max(eqGeometricPSFSlice{end}),'b','LineWidth',geometricSliceFigParams.lineWidth);
plot([eqGeometricBlurCircleDiameterMm(1) eqGeometricBlurCircleDiameterMm(1)],[0 0.5],'r:','LineWidth',geometricSliceFigParams.lineWidth-1);
plot([eqGeometricBlurCircleDiameterMm(end) eqGeometricBlurCircleDiameterMm(end)],[0 0.5],'b:','LineWidth',geometricSliceFigParams.lineWidth-1);
xlabel('Retinal Radius (mm)','FontSize',geometricSliceFigParams.labelFontSize);
ylabel('Point Spread Function','FontSize',geometricSliceFigParams.labelFontSize);
title('Pinhole Camera - Geometric Optics Limited Blur','FontSize',geometricSliceFigParams.titleFontSize);
cbFigAxisSet(geometricSliceFig,geometricSliceFigParams);
legend({sprintf('Pupil: %0.2f mm',pupilDiametersMm(1)) sprintf('Pupil: %0.2f mm',pupilDiametersMm(end))},'Location','NorthEast','FontSize',geometricSliceFigParams.legendFontSize);
FigureSave('PinholeOpticsBlurgeometricSlice',geometricSliceFig,geometricSliceFigParams.figType);







