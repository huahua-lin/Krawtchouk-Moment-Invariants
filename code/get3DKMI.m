function invariants = get3DKMI(arr, dim, varargin)
% get3DKMI Extract the Krawtchouk Moment Invariants (KMI) for 3D images
%==========================================================================
%
% USAGE         invariants = get3DKMI(arr, dim)
%        or...  invariants = get3DKMI(arr, dim, POI)
%        or...  invariants = get3DKMI(arr, dim, order)
%        or...  invariants = get3DKMI(arr, dim, POI, order)
%
% INPUTS
%
%               arr    - Mandatory    - NxNxN array - An array defining the 3D image.
%
%               dim    - Mandatory    - int         - An integer defining the intended size of the arr.
%
%               POI    - Optional     - 1x3 array   - The x,y,z coordinates of point-of-interest.
%
%               order  - Optional     - int         - The order of the moments.
%
% OUTPUTS
%
%               invariants  - 1xP array   - Krawtchouk moment invariants
%==========================================================================
if ndims(arr) ~= 3
    error('Incorrect dimensions. Please input a 3D array.')
end

if dim <= 1
    error('Incorrect value of the size.')
end

if mod(dim,2)==1
    dim = dim-1;  % dim is an even number
end

resizeScale = 1;
if any(size(arr) ~= dim)
    resizeScale = dim / max([size(arr,1), size(arr,2), size(arr,3)]);
    scaledArr = imresize3(arr, resizeScale);
    img = zeros(dim, dim, dim);
    img(1:size(scaledArr,1), 1:size(scaledArr,2), 1:size(scaledArr,3)) = scaledArr;
else
    img = arr;
end

if length(varargin) > 2
    error('Too many inputs.')
end

local = false;
relPOI = [];
order = 2;
for idx = 1:length(varargin)
    if size(varargin{idx}, 2) == 3
        POI = varargin{idx};
        local = true;
        relPOI = POI * resizeScale;
        if any(relPOI < 0) || any(relPOI > size(img, 1))
            error('Incorrect value of the point-of-interest.')
        end

    elseif size(varargin{idx}, 2) == 1
        order = varargin{idx};
        if order < 0
            error('The order must be greater than 0.')
        end
    end
    
end

if local
    disp('Local features are being extracted...')
else
    disp('Global features are being extracted...')
end
    
const = prepStep(dim, order);
invariants = extractFeatures(img, local, relPOI, const);
disp('The extracted features are:')
disp(invariants)
