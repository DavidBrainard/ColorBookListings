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

% Wiki target
wikiTargetDir = '/Users/Shared/GitWebSites/ColorBookListings.wiki';
wikiBaseURL = 'http://davidbrainard.github.io/ColorBookListings/';
wikiTemplateFile = 'WikiHomeTemplate.md';
wikiAddedFile = 'WikiHomeAdded.md';
wikiTargetFile = 'Home.md';
unix(['rm -rf ' fullfile(wikiTargetDir,wikiAddedFile)]);
wikiFid = fopen(fullfile(wikiTargetDir,wikiAddedFile),'w');

% OpticsImage
nSections = nSections + 1;
cbSections{nSections} = 'cbOpticsImage';
cbSectionNames{nSections} = 'Optics and Image';
sectionListings{nSections} = {...
    'cbPinholeOpticsBlur' ...
    };

%% Loop and do each listing in each section
listingDir = pwd;
for s = 1:nSections
    cd(cbSections{s});
    fprintf(wikiFid,['\n## ' cbSectionNames{s} '\n']);
    sectionWebDir = fullfile(webTargetDir,cbSections{s},'');
    if (~exist(sectionWebDir,'dir'))
        mkdir(sectionWebDir);
    end
    sectionDir = pwd;
    for l = 1:length(sectionListings{s})
        cd(sectionListings{s}{l});
        publish(sectionListings{s}{l}); close all;
        sourceDir = 'html';
        targetDir = fullfile(webTargetDir,cbSections{s},sectionListings{s}{l},'');
        unix(['rm -rf ' targetDir]);
        unix(['mv ' sourceDir ' ' targetDir]);
        
        % Add to wiki file
        summaryText = cbGetSummaryText([sectionListings{s}{l} '.m']);
        fprintf(wikiFid,['* [' sectionListings{s}{l} '](' wikiBaseURL cbSections{s} '/' sectionListings{s}{l} '/' sectionListings{s}{l} '.html) - ' summaryText '\n']);
        cd(sectionDir);
    end
    cd(listingDir);
end

%% Close addendum file
fclose(wikiFid);

%% Concatenate to make new wiki page, call git to try to push everything.
% Not sure whether this will handle additions properly.
% 
% And need to figure out how to push without the need for a password
% prompt.  Using the git keychain access seems to work from the command
% line, but not from here.
%   git config --global credential.helper osxkeychain
cd(wikiTargetDir);
unix(['cat ' wikiTemplateFile ' ' wikiAddedFile ' > ' wikiTargetFile]);
unix(['git commit -a -m "Autopublish"; ~/bin/runGitPush.sh']);
cd(webTargetDir);
unix(['git commit -a -m "Autopublish"; ~/bin/runGitPush.sh']);
cd(listingDir);






