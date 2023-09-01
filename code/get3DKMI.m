function invariants = get3DKMI(mask, S, varargin)
% get3DKMI Extract the Krawtchouk Moment Invariants (KMI) for 3D images
%==========================================================================
% AUTHOR        Huahua Lin
% CONTACT       hl4u21@soton.ac.uk
% INSTITUTION   The University of Southampton
% DATE          August 2023
%
% USAGE         invariants = get3DKMI(mask, S)
%        or...  invariants = get3DKMI(mask, S, poi)
%        or...  invariants = get3DKMI(mask, S, order)
%        or...  invariants = get3DKMI(mask, S, poi, order)
%
% INPUTS
%
%               mask   - Mandatory   - NxNxN array  - An array defining the 3D image.
%
%               S      - Mandatory   - int          - An integer defining the intended size of the mask.
%
%               poi    - Optional     - 1x3 array   - The x,y,z coordinates of point-of-interest.
%
%               order  - Optional     - int         - The order of the moments.
%
% OUTPUTS
%
%               invariants  - 1xP array   - Krawtchouk moment invariants
%==========================================================================
if ndims(mask) ~= 3
    error('Incorrect dimensions. Please input a 3D array.')
end

if S <= 1
    error('Incorrect value of the size.')
end

if mod(S,2)==1
    S = S-1;  % S is an even number
end

resizeScale = 1;
if any(size(mask) ~=S )
    resizeScale = S / max([size(mask,1), size(mask,2), size(mask,3)]);
    scaledMask = imresize3(mask, resizeScale);
    img = zeros(S, S, S);
    img(1:size(scaledMask,1), 1:size(scaledMask,2), 1:size(scaledMask,3)) = scaledMask;
else
    img = mask;
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

const = prepStep(S, order);
invariants = extractFeatures(img, local, relPOI, const);
