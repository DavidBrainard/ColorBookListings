% cbPublishListings
%
% Master script to publish all the listings and move the html to
% the website folder.
%
% (c) David Brainard and Andrew Stockman, 2014

%% Clear
clear; close all;

%% Parameters

% Web site target
webTargetDir = '/Users/Shared/GitWebSites/ColorBookListings';
cbSections = {};
nSections = 0;

% OpticsImage
nSections = nSections + 1;
cbSections{nSections} = 'cbOpticsImage';
sectionListings{nSections} = {...
    'cbPinholeOpticsBlur' ...
    };

%% Loop and do each listing in each section
listingDir = pwd;
for s = 1:nSections
    cd(cbSections{s});
    sectionWebDir = fullfile(webTargetDir,cbSections{s},'');
    if (~exist(sectionWebDir,'dir'))
        mkdir(sectionWebDir);
    end
    sectionDir = pwd;
    for l = 1:length(sectionListings{s})
        cd(sectionListings{s}{l});
        publish(sectionListings{s}{l});
        sourceDir = 'html';
        targetDir = fullfile(webTargetDir,cbSections{s},sectionListings{s}{l},'');
        unix(['rm -rf ' targetDir]);
        unix(['mv ' sourceDir ' ' targetDir]);
        cd(sectionDir);
    end
    cd(listingDir);
end

