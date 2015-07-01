function cbAngleToMm
% cbAngleToMm
%
% Illustrates calculations converting degrees of visual angle to mm, by
% making plots of eccentricity in mm versus eccentricity in degrees and the
% other way around.
%
% The underlying routines are in the Psychophysics Toolbox.
%
% Requires: Psychophysics Toolbox
%
% (c) David Brainard and Andrew Stockman, 2015

%% Clear
clear; close all;

%% Set parameters
%
% Can calculate with parameters set either for
% 'Human' or for 'Rhesus' eye.
species = 'Human';
switch (species)
    case 'Human'
        eyeLengthMm = EyeLength(species,'Rodieck');
    case 'Rhesus'
        eyeLengthMm = EyeLength(species,'PerryCowey');              
    otherwise
        error('Unknown species specified');
end

%% Simple trigonometric calculation
degPerMm = DegreesToRetinalMM(1,eyeLengthMm);
mmPerDeg = RetinalMMToDegrees(1,eyeLengthMm);
fprintf('Computing for species %s\n',species);
fprintf('\tMm per degree, linear approximation: %0.2f, assuming nodal point to retina of %0.1f mm\n',degPerMm,eyeLengthMm);
fprintf('\tCorresponding degrees per mm: %0.1f\n',mmPerDeg);
tolerance = 1e-6;
if (abs(degPerMm*mmPerDeg-1) > tolerance)
    error('Linear calculation does not self-invert');
end

%% Degrees to Mm calculations
upperLimitDegrees = 100;
eccDegrees = linspace(0,upperLimitDegrees,100);
eccMm = DegreesToRetinalEccentricityMM(eccDegrees,species,'DaceyPeterson');
eccMmLinear = DegreesToRetinalEccentricityMM(eccDegrees,species,'Linear',eyeLengthMm);
eccDegreesCheck =  RetinalEccentricityMMToDegrees(eccMm,species,'DaceyPeterson');
if (any(abs(eccDegrees-eccDegreesCheck) > tolerance))
    error('Conversion functions do not self invert');
end


%% Mm to degrees calculations
upperLimitMm = 20;
eccMm1 = linspace(0,upperLimitMm,100);
eccDegrees1 = RetinalEccentricityMMToDegrees(eccMm1,species,'DaceyPeterson');
eccDegreesLinear1 = RetinalEccentricityMMToDegrees(eccMm1,species,'Linear',eyeLengthMm);
eccMmCheck1 = DegreesToRetinalEccentricityMM(eccDegrees1,species,'DaceyPeterson');
if (any(abs(eccMm1-eccMmCheck1) > tolerance))
    error('Conversion functions do not self invert');
end

% This is the tangent calculation, which although it once seemed to me
% a better approximation than the pure linear approximation turns out to
% deviate more from the model eye based calculation.  Not plotted but you
% could be commenting out the corresponding lines in the plot code below.
eccMmTangent = DegreesToRetinalMM(eccDegrees,eyeLengthMm,true);
eccDegreesTangent1 = RetinalMMToDegrees(eccMm1,eyeLengthMm,true);

%% Drasdo and Fowler Figure 2
%
% This is ellipsoid model data, digitized from their figure.  It is not
% clear to me which curve (ellipsoid or sphere) Dacey and Peterson
% digitized and fit.  First column is degrees, second is mm.
drasdoFowlerData = ...
    [1.429393919263E0	4.001594366571E-1
    3.492105671073E0	9.230850993158E-1
    6.034234593326E0	1.662020859629E0
    1.000354303683E1	2.739148342523E0
    1.571668991785E1	4.215850660998E0
    2.381363626298E1	6.400797183286E0
    3.286388095396E1	8.862844615691E0
    4.016364401337E1	1.073941407028E1
    4.572178303328E1	1.227837640337E1
    4.794503864125E1	1.289396133661E1
    4.984831373591E1	1.335521158573E1
    5.254323612126E1	1.396991961735E1
    5.523705130760E1	1.455364379194E1
    6.126021391085E1	1.590582608118E1
    6.632675657123E1	1.688795588919E1
    7.012555636750E1	1.759356938816E1
    7.503155517173E1	1.848303992560E1
    8.135698310414E1	1.949381518634E1
    8.689076374588E1	2.035113266459E1
    9.321065568325E1	2.120698864014E1
    9.921278150535E1	2.197047764565E1];

%% Make a figure for the box
%
% By uncommenting the plot line for the Drasdo and Fowler data
% you can see that the approximation provided by the fit is so-so,
% although the formula is much closer to the Drasdo and Fowler curve 
% than the linear approximation.  I don't think the juice is worth the
% squeeze as far as trying to do better with that figure, given individual
% differences and model eye differences.  If one wanted to go down that
% road, probably starting with a ray trace of our best current model eye
% would be the place to start.
[angleToMmFig,figParams] = cbFigInit;
figParams.xLimLow = 0;
figParams.xLimHigh = 80;
figParams.xTicks = [0 10 20 30 40 50 60 70 80];
figParams.xTickLabels = {};
figParams.yLimLow = 0;
figParams.yLimHigh = 25;
figParams.yTicks = [0 5 10 15 20 25];
figParams.yTickLabels = {};

plot(eccDegrees,eccMmLinear,'b:','LineWidth',figParams.lineWidth);
plot(eccDegrees,eccMm,'r','LineWidth',figParams.lineWidth);
% plot(eccDegrees,eccMmTangent,'g','LineWidth',figParams.lineWidth);
plot(drasdoFowlerData(:,1),drasdoFowlerData(:,2),'k:','LineWidth',figParams.lineWidth-1);

xlabel('Eccentricity (degrees)','FontSize',figParams.labelFontSize);
ylabel('Eccentricity (mm)','FontSize',figParams.labelFontSize);
title('Eccentricity Conversion','FontSize',figParams.titleFontSize);
cbFigAxisSet(angleToMmFig,figParams);
legend({'Linear', 'Model Eye Based'},'Location','NorthWest','FontSize',figParams.legendFontSize);
FigureSave('EccentricityAngleToMmConversion',angleToMmFig,figParams.figType);

%% Can also make the inverse figure.
[mmToAngleFig,figParams] = cbFigInit;
figParams.xLimLow = 0;
figParams.xLimHigh = 25;
figParams.xTicks = [0 5 10 15 20 25];
figParams.xTickLabels = {};
figParams.yLimLow = 0;
figParams.yLimHigh = 80;
figParams.yTicks = [0 10 20 30 40 50 60 70 80];
figParams.yTickLabels = {};

plot(eccMm1,eccDegreesLinear1,'b:','LineWidth',figParams.lineWidth);
plot(eccMm1,eccDegrees1,'r','LineWidth',figParams.lineWidth);
% plot(eccMm1,eccDegreesTangent1,'g','LineWidth',figParams.lineWidth);

xlabel('Eccentricity (mm)','FontSize',figParams.labelFontSize);
ylabel('Eccentricity (degrees)','FontSize',figParams.labelFontSize);
title('Eccentricity Conversion','FontSize',figParams.titleFontSize);
cbFigAxisSet(mmToAngleFig,figParams);
legend({'Linear', 'Model Eye Based'},'Location','NorthWest','FontSize',figParams.legendFontSize);
FigureSave('EccentricityMmToAngleConversion',mmToAngleFig,figParams.figType);


