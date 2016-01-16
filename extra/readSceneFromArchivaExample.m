

%% Load in an image
%
% Grab a hyperspectral scene from our archiva server
% and go from there.

% Load a scene and extract a normalized grayscale image
client = RdtClient('isetbio');
remotePath = 'resources/scenes/hyperspectral/manchester_database';
client.crp(remotePath);
[theScene, theSceneArtifact] = client.readArtifact('scene7', 'type', 'mat');
imageData = theScene.scene.data.photons(:,:,15);

% Resize
data.nImagePixels = 512;
imageData = imresize(imageData,[data.nImagePixels,data.nImagePixels]);

% Quick and dirty tone map.  The image looked a little dark. 
%   Scale down, take square root, normalize, and square 
%   makes it look a bit better.
downFactor = 1.5;
imageData = imageData/(downFactor*max(imageData(:)));
imageData = sqrt(imageData);
imageDdata = imageData.^2;
data.imageData = imageData/max(imageData(:));
imshow(data.imageData);