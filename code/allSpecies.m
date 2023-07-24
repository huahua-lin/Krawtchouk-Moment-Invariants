S = 128; order = 2;
data = readtable('../samples/2023-04-25_Species_Classifications.xlsx');
numOfData = size(data,1);
disp(numOfData);
Desc = zeros(numOfData,3);
for i = 1:numOfData
    tic;
    disp(i);
    dataLocation = string(data{i,5});
    folder = '../CP40_AI_50/' + dataLocation(1,:) + '_recon/';
    FileList = natsortfiles(dir(fullfile(folder, '*.png')));
    rows = 512; cols = 512;
    numOfSlices = length(FileList);
    
    array3d = zeros(rows, cols, numOfSlices);
    for slice = 1 : numOfSlices
        File = fullfile(FileList(slice).folder, FileList(slice).name);
        thisSlice = imread(File);
        array3d(:,:,slice) = thisSlice;
    end
    
    mask = false(size(array3d));
    pixIndex = find(array3d > 0);
    mask(pixIndex) = 1;
    
    resizeScale = 128 / max([rows, cols, numOfSlices]);
    scaledMask = imresize3(mask, resizeScale);

    gridOri = zeros(S, S, S);
    gridOri(1:size(scaledMask,1), 1:size(scaledMask,2), 1:size(scaledMask,3)) = scaledMask;
    
    const = prepStep(S,order);
    desc1 = extractFeatures(gridOri, const);

    Desc(i,:) = desc1;
end

