% Set defaults for EphysViewer
% rangequest=1 lets you decide what range of data you want to import
% qmquestion=1 lets you import the data in a slower but more memory efficient manner
% chunkmode=1 lets you split the data into chunks based on some trigger
% channels=1 lets you chose only certain channels to import
% filter=x runs medfilt1 on the data with order x

% Adam Packer
% April 10th, 2007
% Added filter capability December 18th, 2009


function [rangequest,qmquestion,chunkmode,channels,filter]=ephysdefaults

rangequest=0;
qmquestion=0;
chunkmode=1;
channels=0;
filter=0;