function cbPinholeOpticsBlur
% cbAngleToMm
%
% Illustrates calculations converting degrees of visual angle to mm.
%
% Requires: Psychophysics Toolbox
%
% (c) David Brainard and Andrew Stockman, 2015

%% Clear
clear; close all;

%% Set parameters
eyeNodalToRetinaMm = 17;

%% Simple trigonometric calculation
degPerMm = eyeNodalToRetinaMm*tan(deg2rad(1));
mmPerDeg = 1/degPerMm;
fprintf('Mm per degree: %0.2f, assuming nodal point to retina of %0.1f mm\n',degPerMm,eyeNodalToRetinaMm);
fprintf('\tCorresponding degrees per mm: %0.1f\n',mmPerDeg);

%% You can get fancier, by ray tracing a model eye and then approximating what you find.
%
% Drasdo and Fowler, 1974 (British J. Opthth, 58,pp. 709 ff. report the
% results of such an enterprise.  They give the formula
%
%   T = 11.06*sin(5.181*d)/[53.7*sin(theta)]
%
% where d is the polar eccentricity of the point in mm of retina being considerd and T is the
% number of tangential mm correspdonding to 1 degree at that eccentricity.
%
% This formula allows us to work out the number of mm per degree as a
% function of mm, and we can integrate this to get a curve which plots mm
% versus degrees.

dStepSizeMm = 0.1;
dStartMm = 0;

