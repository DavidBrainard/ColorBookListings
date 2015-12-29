function varargout = cbOpticsImage_HumanPSF(varargin)
%
% Illustrate human point spread functions.
%
% Uses population wavefront aberration data to make PSF images
% for a series of pupil sizes.
%
% (c) David Brainard and Andrew Stockman, 2015

varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams) 

%% Hello
clear global; ieInit;
UnitTest.validationRecord('SIMPLE_MESSAGE', sprintf('%s',mfilename));
outputDir = sprintf('%s_Output',mfilename);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
close all; drawnow;

%% The tutorial only uses 1 wavelength at a time. So, for plotting, we use
% this index.
maxMM = 2;
maxUM = 20;      
pupilfuncrangeMM = 5;

%% Start by making pictures for a diffraction limited PSF
%
% Use wvfCreate to create a wavefront variable to explore with.
%
% This wavefront by default has the 0's for all zernike coeffs
% Notice that the calcpupilMM is by default 3, meaning we are simulating
% the wavefront PSF for a pupil of 3mm diameter.
wvf0 = wvfCreate;

% Look at the plot of the normalized PSF within 1 mm of the center.
% Variable maxUM is used to specify size of plot from center of the PSF.
%
% The plot shows an airy disk computed from the Zernike polynomials; that
% is representing the diffraction-limited PSF obtained when the Zernike
% coefficients are all zero.
maxMinutes = 2;
wl = wvfGet(wvf0,'calc wavelengths');

calcp = 2;
wvf0 = wvfSet(wvf0,'calc pupil size',calcp);
wvf0 = wvfComputePSF(wvf0);
figure;
subplot(4,3,1);
wvfPlot(wvf0,'2d psf angle','min',wl,maxMinutes,'no window');
title(sprintf('%d nm, %d mm pupil',wl,calcp));

calcp = 4;
wvf0 = wvfSet(wvf0,'calc pupil size',calcp);
wvf0 = wvfComputePSF(wvf0);
subplot(4,3,4);
wvfPlot(wvf0,'2d psf angle','min',wl,maxMinutes,'no window');
title(sprintf('%d nm, %d mm pupil',wl,calcp));

calcp = 6;
wvf0 = wvfSet(wvf0,'calc pupil size',calcp);
wvf0 = wvfComputePSF(wvf0);
subplot(4,3,7);
wvfPlot(wvf0,'2d psf angle','min',wl,maxMinutes,'no window');
title(sprintf('%d nm, %d mm pupil',wl,calcp));

calcp = 8;
wvf0 = wvfSet(wvf0,'calc pupil size',calcp);
wvf0 = wvfComputePSF(wvf0);
subplot(4,3,10);
wvfPlot(wvf0,'2d psf angle','min',wl,maxMinutes,'no window');
title(sprintf('%d nm, %d mm pupil',wl,calcp));

%% Wavefront measurements of human eyes and the effects of single-vision
% corrective eyeglasses 
%
% We have access to measurements of the pupil function of real human eyes. The
% optics of these eyes are not perfect, so they have interesting pupil functions
% and PSF shapes.

% Set up the wvf structure
measMM = 6;
calcMM = 3;
maxMM = 3;
theWavelengthNM = 550;
wvfHuman0 = wvfCreate('measured pupil',measMM,'calculated pupil',calcMM);
wvfHuman0 = wvfSet(wvfHuman0,'wavelength',theWavelengthNM);

% Load in some measured data
sDataFile = fullfile(wvfRootPath,'data','sampleZernikeCoeffs.txt');
theZernikeCoeffs = importdata(sDataFile);
whichSubjects = [3 7];
theZernikeCoeffs = theZernikeCoeffs(:,whichSubjects);
nSubjects = size(theZernikeCoeffs,2);
nRows = ceil(sqrt(nSubjects));
nCols = ceil(nSubjects/nRows);

% Plot subject PSFs, one by one
for ii = 1:nSubjects
    fprintf('** Subject %d\n',whichSubjects(ii))

    wvfHuman = wvfSet(wvfHuman0,'zcoeffs',theZernikeCoeffs(:,ii));
    wvfHuman = wvfComputePSF(wvfHuman);
    
    vcNewGraphWin;
    subplot(2,2,1);
    wvfPlot(wvfHuman,'2d pupil amplitude space','mm',[],calcMM,'no window');
    subplot(2,2,2);
    wvfPlot(wvfHuman,'2d pupil phase space','mm',[],calcMM,'no window');
    subplot(2,2,3:4);
    wvfPlot(wvfHuman,'2d psf space','mm',[],maxMM,'no window');
end

%% Single-vision eyewear generally corrects only the lowest-order
% Zernike aberrations (defocus given in diopters) and astigmatism (cylinder
% correction also given in diopters). The Zernike coefficients give us an
% easy and convenient way to simulate corrective lenses; we can simply set
% those Zernike coefficients to zero and see what the PSFs look like!
%
% Plot their corrected PSFs, one by one, How do the corrected PSFs compare
% to the uncorrected ones? their peaks? their widths?
%
% Try changing the whichSubjects array above to look at other sample data. Do
% eyeglasses help correct the aberrations in those subjects?
%
% If you were to spend thousands of dollars on laser eye surgery, would you
% want them to only correct the first order of wavefront aberrations, like
% eyeglasses, or do a full wavefront measurement?
% 
% Suppose you knew that such surgery would correct some of the lower order aberrations but some of the
% higher order aberrations worse.  How would you compute the net effect of
% something like that?
for ii = 1:nSubjects
    fprintf('** Subject %d corrected\n',whichSubjects(ii))
    
    % Correct defocus and astigmatism
    zCoeffs = theZernikeCoeffs(:,ii);
    zCoeffs(4:6) = 0;
    wvfHuman = wvfSet(wvfHuman0,'zcoeffs',zCoeffs);
    wvfHuman = wvfComputePSF(wvfHuman);
    
    vcNewGraphWin;
    subplot(2,2,1);
    wvfPlot(wvfHuman,'2d pupil amplitude space','mm',[],calcMM,'no window');
    subplot(2,2,2);
    wvfPlot(wvfHuman,'2d pupil phase space','mm',[],calcMM,'no window');
    subplot(2,2,3:4);
    wvfPlot(wvfHuman,'2d psf space','mm',[],maxMM,'no window');
end

end

