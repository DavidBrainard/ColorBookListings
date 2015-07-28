function varargout = cbOpticsImage_PoissonNoiseImages(varargin)
%
% Illustrate magnitude of Poisson noise as a function of light level
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

%% Frozen noise, so that we can validate OK
randomSeedValue = 26;
rng(randomSeedValue);

%% Load a hyperspectral scene in ISETBIO format.
theHyperSceneDir = '/Volumes/Users1/Shared/Matlab/Analysis/hyperspectral-images/manchester_database';
theHyperSceneName = 'isetbioSceneFor_scene7.mat';
theData = load(fullfile(theHyperSceneDir,theHyperSceneName));
vcAddAndSelectObject(theData.scene); sceneWindow;
theData.sceneRGBImage = sceneGet(theData.scene,'rgb image');

%% The scene object has a view about what the illuminant is.
%
% Get this to help us set a reasonable physical scale for the scene.
sceneWave = sceneGet(theData.scene,'wave');
sceneS = WlsToS(sceneWave);
sceneIlluminantSpd = sceneGet(theData.scene,'illuminant energy');
sceneIlluminantXYZ = sceneGet(theData.scene,'illuminant XYZ');
load T_xyz1931;
T_xyz = SplineCmf(S_xyz1931,T_xyz1931,sceneS);
ourSceneIlluminantXYZ = 683*T_xyz*sceneIlluminantSpd*sceneS(2);
UnitTest.assertIsZero(max(abs(sceneIlluminantXYZ(:)-ourSceneIlluminantXYZ(:))),'Check on XYZ computation',0);
fprintf('Scene illumination luminance taken as %0.0f cd/m2\n',sceneIlluminantXYZ(2));

%% Make optical image
%
% We do this to get retinal irradiance, and for the first part of this
% script we only want photon noise so we skip the blurring.
theData.oi = oiCreate('human');
optics = opticsSet(optics,'off axis method','skip');
optics = opticsSet(optics,'otf method','skip otf');
theData.oi = oiCompute(theData.oi,theData.scene);
theData.oiRGBImage = oiGet(theData.oi,'rgb image');
vcAddAndSelectObject(theData.oi); oiWindow;

%% Get the retinal irradiance per pixel
oiWave = oiGet(theData.oi,'wave');
oiIrradiance_PhotonsPerSecM2 = oiGet(theData.oi,'photons');
oiPixelSize_M = oiGet(theData.oi,'sample spacing');
oiPixelArea_M2 = oiPixelSize_M(1)*oiPixelSize_M(2);
oiIrradiance_PhotonsPerSecPixel = oiIrradiance_PhotonsPerSecM2*oiPixelArea_M2;

%% Pick an integration time to get photons
theData.integrationTime_Sec = 0.1;
oiEnergy_PhotonsPerPixel = oiIrradiance_PhotonsPerSecPixel*theData.integrationTime_Sec;

%% Get one wavelength plane
theData.whichWavelength = 550;
wlIndex = find(oiWave == theData.whichWavelength);

%% Make image noise free
meanPhotons0 = double(oiEnergy_PhotonsPerPixel(:,:,wlIndex));
meanScalar = 0.2;
displayGamma = 0.5;
scaledMeanPhotons0 = meanScalar*meanPhotons0/mean(meanPhotons0(:));
figure; clf;
imshow(scaledMeanPhotons0.^displayGamma);
title('No noise')

%% Make images for different irradiance levels
theIrradianceScaleFactors = [1e-6 1e-4 1e-2 1e-1 1];
for ii = 1:length(theIrradianceScaleFactors)
    % Get mean number of photons at each pixel at the desired wavelength
    % and irradiance scale factor
    meanPhotons = theIrradianceScaleFactors(ii)*double(oiEnergy_PhotonsPerPixel(:,:,wlIndex));
    
    % Get a draw of Poisson noise around the mean
    noisyPhotons = poissrnd(meanPhotons);
    
    % Scale in a consistent manner for display
    scaledNoisyPhotons = meanScalar*noisyPhotons/mean(meanPhotons(:));
    
    % Make a figure
    figure; clf;
    imshow(scaledNoisyPhotons.^displayGamma);
    title(sprintf('Poisson photon noise for irradiance factor %g',theIrradianceScaleFactors(ii)));
    
    % Check.  Isetbio can give us a photon noised image, so when the scale
    % factor is unity we generate that too and have a look
    % I AM WORKING HERE.
    % if (theIrradianceScaleF
end


%% Can also make the inverse figure.
if (runTimeParams.generatePlots)
  
end

%% Save validation data
UnitTest.validationData('theData', data);

end


