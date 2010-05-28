% traces = pi_load_data(inputDirectory)
%
% read in tif stacks sequentially. the traces from each neuron will be
% extracted and concatenated across stacks.
%
% each image stack should contain data in the form P x N x T
% where P = pixels per target
%       N = number of targets (putative neurons)
%       T = number of frames per segment
%
% files should be numbered in some reasonable way (eg: 1.tif - 10.tif)
% and do not all need to be the same frame length
%
% the stacks used by this function should probably be generated using the
% compress-stack imagej plugin
% 
% tamachado, Columbia University, 5/18/2010
function traces = pi_load_data(inputDirectory)
currDirectory = pwd;
disp('loading data...')
try
    % get file names
    cd(inputDirectory);
    files = ls;

    % put into a 1D array if needed
    if (min(size(files) > 1))
        f1 = [];
        for ii = 1:size(files,1)
            f1 = [f1 ' ' files(ii,:)]; %#ok<AGROW>
        end
        files = f1;
    end
    
    % get indices of names and spaces
    names = regexp(files,'\d+.tif');
    spaces = regexp(files,'\s');
    off = -1;
    if isequal(computer,'PCWIN') || isequal(computer,'PCWIN64')
        names = regexp(files,'\d+.tif');
        spaces = regexp(files,'\s+');
        files(end+1) = ' ';
        spaces(end+1) = length(files);
        off = 0;
    end
        
    fn = cell(length(names),1);
    for ii = 1:length(names)
        nextSpace = spaces(find(spaces > names(ii),1));
        fn{ii} = files(names(ii):nextSpace+off);
    end
    
    % open each stack and get dimension
    nFrames = zeros(length(names),1);
    for ii = 1:length(names)
        disp(ii)
        inf = imfinfo(fn{ii});
        nFrames(ii) = length(inf);
    end
    nPixels = inf(1).Height;
    nCells  = inf(1).Width;
    traces = zeros(nCells,sum(nFrames));
    
    % process each movie
    for ii = 1:length(names)
        disp(ii)
        % process each frame
        for jj = 1:nFrames(ii)
            traces(:,jj+sum(nFrames(1:ii-1))) = sum(imread(fn{ii},jj),1);
        end
    end
catch err
    disp(err)
    cd(currDirectory);
end
cd(currDirectory);
disp('done!')