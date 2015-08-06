function varargout = cbColorMatching_KonigFundamentals(varargin)
%
% Illustrate the ideas underlying Konig fundamentals.
%
% Shows how to use dichromatic confusion data to estimate cone fundamentals
% from color matching functions.  This is done Stiles-Burch 10 degree cmfs and the
% Stockman-Sharpe 10 degree cone fundamentals, but the principles
% would apply to any tristimulus system and cone fundamentals that
% were a linear transfomration of the color matching functions.
%
% This routine isn't based on real data, the dichromatic confusion data are
% syntehsized using the 'known' cone fundamentals.  What this routine shows
% is how such data lock donw the desired transformation.
%
% The Stiles-Burch 10-degree cmfs are expressed with respect to primaries at
% 645.16, 526.32, 444.44 nm.
%
% See also cbColorMatching_StilesBurch10Cmfs.
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
close all; drawnow;

%% Get and plot Stiles-Burch 10-degree color matching functions
%
% Also spline to 1 nm and 10 nm sampling for plotting purposes.
load T_stiles10;
data.wls = SToWls(S_stiles10);
data.S_stiles10 = S_stiles10;
data.T_stiles10 = T_stiles10;
clear S_stiles10 T_stiles10

% These splines go from 390 to 750, which seems sufficient
S_1nm = [390 1 361];
data.wls_1nm = SToWls(S_1nm);
data.T_stiles10_1nm = SplineCmf(data.wls,data.T_stiles10,S_1nm);
S_10nm = [390 10 37];
data.wls_10nm = SToWls(S_10nm);
data.T_stiles10_10nm = SplineCmf(data.wls,data.T_stiles10,S_10nm);

%% Get the cmf spectrum locus normalized to simplex
for i = 1:size(data.T_stiles10_1nm,2);
    data.T_stiles10_1nm_simplex(:,i) = data.T_stiles10_1nm(:,i)/sum(data.T_stiles10_1nm(:,i));
end
for i = 1:size(data.T_stiles10_10nm,2);
    data.T_stiles10_10nm_simplex(:,i) = data.T_stiles10_10nm(:,i)/sum(data.T_stiles10_10nm(:,i));
end

%% Create the color matching primaries
% Round to nearest nm which is not exact but appears to
% work to within about 1 percent numerically in various
% checks below.
data.B_1nm = zeros(S_1nm(3),3);
data.primaryWls = [645, 526, 444];
for i = 1:3
    wlIndex = find(data.wls_1nm == data.primaryWls(i));
    data.B_1nm(wlIndex,i) = 1;
end

%% Load Stockman-Sharpe 10-degree cone fundamentals
load T_cones_ss10
data.T_cones10_1nm = SplineCmf(S_cones_ss10,T_cones_ss10,data.wls_1nm);

% Fit with linear transform of cmf's, just to show that it works.
data.M_CmfToCones = ((data.T_stiles10_1nm')\(data.T_cones10_1nm'))';
data.T_cones10_fit_1nm = data.M_CmfToCones*data.T_stiles10_1nm;

%% Find the stimuli that isolate each of the cones.
%
% Get the isolating directions.
% We have the transformation for cmfs to cones spectral sensitivies.
% This is also the transformation between tristimulus coordinates
% and cone excitations. Invert this to get transformation between
% cone excitations and tristimulus coordinates.  Then apply to the
% unit cone excitation vectors to get the cone isolating tristimulus
% vectors.
%
% Make a normalized version of the vectors for plotting.
data.M_ConesToCmf = inv(data.M_CmfToCones);
data.coneIsolatingRGBDirs = data.M_ConesToCmf*[[1 0 0]', [0 1 0]', [0 0 1]'];
for i = 1:size(data.coneIsolatingRGBDirs,2)
    data.coneIsolatingRGBDirsNorm(:,i) = data.coneIsolatingRGBDirs(:,i)/norm(data.coneIsolatingRGBDirs(:,i));
end

%% Find the cone sensitivity vectors in RGB tristimulus space.
% These are unit vectors such that when you project onto them, you
% get the L, M, and S responses.  When we work in cone space, thesr
% are just the unit vectors, which are conceptually best expressed
% as row vectors, since then multiplying tristimulus coordinates
% from the left with such a vector directly gives cone excitation.
%
% Since cone excitations themselves must be independent of the space
% we compute them in, we seek a row vector r such that r*t = u*c,
% where r is the response row vector in RGB, t is the tristimulus
% column vector, u is a unit row vector and c is a cone excitation
% column vector.  We also have t = M*c with M being the M_ConesToCmf
% computed above.  Thus we must have r = u*M_inv = u*M_CmfToCones.
% Note that r does not necessarily have unit length -- it is only
% in the cone excitation space where the response vectors are guaranteed
% to be the unit vectors and to have unit length.
%
% Also produce normalized version for plotting.
data.coneResponseRGBVectors = [ [1 0 0] ; [0 1 0] ; [0 0 1] ]*data.M_CmfToCones;
for i = 1:size(data.coneResponseRGBVectors,1)
    data.coneResponseRGBVectorsNorm(i,:) = data.coneResponseRGBVectors(i,:)/norm(data.coneResponseRGBVectors(i,:));
end

%% Get the cmf spectrum locus normalized to simplex (R + G + B = 1)
for i = 1:size(data.T_stiles10_1nm,2);
    data.T_stiles10_1nm_simplex(:,i) = data.T_stiles10_1nm(:,i)/sum(data.T_stiles10_1nm(:,i));
end
for i = 1:size(data.T_stiles10_10nm,2);
    data.T_stiles10_10nm_simplex(:,i) = data.T_stiles10_10nm(:,i)/sum(data.T_stiles10_10nm(:,i));
end

%% Get the cone isolating dirs normalized to the simplex
for i = 1:size(data.coneIsolatingRGBDirs,2);
    data.coneIsolatingRGBDirs_simplex(:,i) = data.coneIsolatingRGBDirs(:,i)/sum(data.coneIsolatingRGBDirs(:,i));
end

%% Get dichromatic confusion lines
%
% We can generate these by adding the cone isolating direction to
% to any stimulus.  So let's add it to each stimulus on the spectrum
% locus.
%
% These converge on the chromaticity of the isolating direction for each
% My intuition for this is that as you add more and more of the stimulus in
% the cone isolating direction, it dominates the tristimulus coordinates
% more and more, swamping whatever it was being added to.  In the limit,
% then, the chromaticity of the summed stimulus will be that of the cone
% isolating stimulus.
%
% To make these plot nicely, we chose a scale factor for each type of
% dichromat differently.  For the M cone isolating direction, it's the
% negative excursion that intersects the simplex; plotting the confusion
% lines looks nicer if we use the negative rather than the positive
% excusion.
for w = 1:3
    switch (w)
        % Protanope
        case 1
            whichConfusionColor = 'r';
            confusionLineLengthFactor = 1;
            % Deuteranope
        case 2
            whichConfusionColor = 'g';
            confusionLineLengthFactor = -3;
            % Tritanope
        case 3
            whichConfusionColor = 'b';
            confusionLineLengthFactor = 30;
    end
    for i = 1:size(data.T_stiles10_10nm,2);
        nConfusionPoints = 100;
        for j = 1:nConfusionPoints
            confusionLine{w,i}(:,j) = data.T_stiles10_10nm(:,i) + ...
                confusionLineLengthFactor*((j-1)/nConfusionPoints)*data.coneIsolatingRGBDirs(:,w);
            confusionLine_simplex{w,i}(:,j) = confusionLine{w,i}(:,j)/sum(confusionLine{w,i}(:,j));
        end
    end
end

%% If we have measured the confusion lines, we know the cone isolating directions.
%
% This is trivial in the case that we measure them in the full tristimulus
% space -- we just need to find the direction of the confusion line for any
% base stimulus and we have it.  As we'll show below, having the direction
% of the cone isolating stimulus for each cone class is enough to lock
% down the transformation from tristimulus coordinates to cone exciations,
% and this together with the color matching functions is enough to give us
% the cone fundamentals.  There is a free scaling parameter left for each
% fundamental, that has to be locked down some other way.
%
% What if we have just the chromaticities of several confusion lines for
% each type of dichromat.  Then we can still do what we need.  This isn't
% too hard, but is a little less trivial.  We use the confusion lines for
% each type of dichromat to find where they intersect.  This "copunctal point"
% gives us the chromaticity of the cone isolating directions, and from there
% we can get the tristimulus coordinates the cone isolating directions up
% to a free scale factor.  Let's illustrate that.
%
% We can express each confusion line as g = a*r+b.   First step is to find
% a and b for each confusion line.

%% Plot spectrum locus and isolating vectors in the r-g chromaticity plane
%
% This includes a little empirical calculation of the gamut we can obtain
% with positive combinations of the isolating vectors.  This is a bit
% non-intuitive, to me at least.
[chromaticityFig,figParams] = cbFigInit;
set(gcf,'Position',[100 254 1163 446]);
if (runTimeParams.generatePlots)
    figParams.xLimLow = -2.5;
    figParams.xLimHigh = 1.5;
    figParams.xTicks = [-2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5];
    figParams.xTickLabels = {'^{ }-2.5_{ }' '^{ }-2.0_{ }' '^{ }-1.5_{ }' '^{ }-1.0_{ }' '^{ }-0.5_{ }' '^{ }0.0_{ }' ...
        '^{ }0.5_{ }' '^{ }1.0_{ }' '^{ }1.5_{ }'};
    figParams.yLimLow = -1;
    figParams.yLimHigh = 3;
    figParams.yTicks = [-1 -0.5 0 0.5 1 1.5 2 2.5 3.0];
    figParams.yTickLabels = {'^{ }-1.0_{ }' '^{ }-0.5_{ }' '^{ }0.0_{ }' ...
        '^{ }0.5_{ }' '^{ }1.0_{ }' '^{ }1.5_{ }' '^{ }2.0_{ }' '^{ }2.5_{ }' '^{ }3.0_{ }'};
    
    % Montage plot of confusion lines
    for w = 1:3
        subplotHandle = subplot(1,3,w); hold on;
        set(gca,'FontName',figParams.fontName,'FontSize', ...
            figParams.axisFontSize-figParams.subplotFontShrink,'LineWidth',figParams.axisLineWidth);
        
        % Plot the confusion lines
        switch (w)
            case 1
                % Protanope
                whichConfusionColor = 'r';
                titleStr = 'Protan';
            case 2
                % Deuteranope
                whichConfusionColor = 'g';
                titleStr = 'Deutan';
            case 3
                % Tritanope
                titleStr = 'Tritan';
                whichConfusionColor = 'b';
        end
        for i = 1:size(data.T_stiles10_10nm,2);
            plot(confusionLine_simplex{w,i}(1,:),confusionLine_simplex{w,i}(2,:),whichConfusionColor,'LineWidth',1);
        end
        
        % Plot the spectrum locus on the diagram along with equal energy white.
        plot(data.T_stiles10_1nm_simplex(1,:)',data.T_stiles10_1nm_simplex(2,:)', ...
            'k','LineWidth',figParams.lineWidth);
        plot(data.T_stiles10_10nm_simplex(1,:)',data.T_stiles10_10nm_simplex(2,:)', ...
            'ko','MarkerFaceColor','k','MarkerSize',figParams.markerSize-14);
        
        % Plot where the cone isolating dirs lie on the diagram
        %
        % The M-cone chromaticity corresponds to the negative direction
        % of the primary, so it's plotted without a fill.
        plot([data.coneIsolatingRGBDirs_simplex(1,1)], ...
            [data.coneIsolatingRGBDirs_simplex(2,1)], ...
            'ro','MarkerFaceColor','r','MarkerSize',figParams.markerSize-10);
        plot([data.coneIsolatingRGBDirs_simplex(1,2)], ...
            [data.coneIsolatingRGBDirs_simplex(2,2)], ...
            'go','MarkerSize',figParams.markerSize-10);
        plot([data.coneIsolatingRGBDirs_simplex(1,3)], ...
            [data.coneIsolatingRGBDirs_simplex(2,3)], ...
            'bo','MarkerFaceColor','b','MarkerSize',figParams.markerSize-10);
        
        % Labels
        xlabel('r','FontSize',figParams.labelFontSize-figParams.subplotFontShrink);
        ylabel('g','FontSize',figParams.labelFontSize-figParams.subplotFontShrink);
        title(titleStr,'FontSize',figParams.titleFontSize-figParams.subplotFontShrink);
        axis('square');
        cbFigAxisSet(subplotHandle,figParams);
    end
    
    % Save the figure
    FigureSave(fullfile(outputDir,[mfilename '_ConfusionLines_rgChrom']),chromaticityFig,figParams.figType);
end

%% Save validation data
UnitTest.validationData('validateDataStruct', data);

end


