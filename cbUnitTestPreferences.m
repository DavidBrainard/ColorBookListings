% Method to set cb-specific unit testing preferences. 
%
% Generally, this function should be edited for your site and then run
% once.
%
% To use the unit testing framework in different projects, copy this file
% to the projects' validation directory and adapt the p-struct according
% to that project's specifics.  You will want to change isetbioValidation
% and isetbioRootPath and so forth.
%
% If you just run the distributed (Template) version to accept the
% defaults, the FAST validation in isetbio will run.  But the FULL forms of
% validation need the validation data. These can be found here:
%
%  http://????? - To be determined.
%
%
% ISETBIO Team, 2015

function cbUnitTestPreferencesTemplate

    % Specify project-specific preferences
    p = struct(...
            'projectName',           'cbScripts', ...                                                                                 % The project's name (also the preferences group name)
            'validationRootDir',     '/Users/Shared/Matlab/Analysis/ColorBookListings', ...                                           % Directory location where the 'scripts' subdirectory resides.
            'alternateFastDataDir',  '',  ...                                                                                         % Alternate FAST (hash) data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/fast
            'alternateFullDataDir',  '/Volumes/Users1/Shared/Matlab/Analysis/ColorBookListings', ...                                  % Alternate FULL data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/full
            'clonedWikiLocation',    '/Users/Shared/GitWebSites/ColorBookListings.wiki', ...                                          % Local path to the directory where the wiki is cloned. Only relevant for publishing tutorials.
            'clonedGhPagesLocation', '/Users/Shared/GitWebSites/ColorBookListings', ...                                               % Local path to the directory where the gh-pages repository is cloned. Only relevant for publishing tutorials.
            'githubRepoURL',         'https://github.com/DavidBrainard/ColorBookListings', ...                                        % Github URL for the project. This is only used for publishing tutorials.
            'generateGroundTruthDataIfNotFound',   true ...                                                                           % Flag indicating whether to generate ground truth if one is not found
        );

    generatePreferenceGroup(p);
    UnitTest.usePreferencesForProject(p.projectName);
    
end

function generatePreferenceGroup(p)
    % remove any existing preferences for this project
    if ispref(p.projectName)
        rmpref(p.projectName);
    end
    
    % generate and save the project-specific preferences
    setpref(p.projectName, 'projectSpecificPreferences', p);
    fprintf('Generated and saved preferences specific to the ''%s'' project.\n', p.projectName);
end