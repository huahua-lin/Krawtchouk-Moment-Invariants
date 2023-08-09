function invariants = KMI3D(img, varargin)
% KMI3D Extrat the Krawtchouk Moment Invariants (KMI) for 3D images
%==========================================================================
% AUTHOR        Huahua Lin
% CONTACT       hl4u21@soton.ac.uk
% INSTITUTION   The University of Southampton
% DATE          August 2023
%
% USAGE         features = 3DKMI(img)
%        or...  features = 3DKMI(img, poi)
%        or...  features = 3DKMI(img, order)
%        or...  features = 3DKMI(img, poi, order)
% INPUTS        
%
%               img    - Mandatory    - NxNxN array - An array defining the 3D image.
%
%               poi    - Optional     - 1x3 array   - The x,y,z coordinates of point-of-interest.                        
%
%               order  - Optional     - int         - The order of the moments.
%
% OUTPUTS
%
%               invariants  - 1xP array   - Krawtchouk moment invariants
%==========================================================================

poi = zeros(1,3);
order = 0;

if size(varargin) == 0
    if ndims(img) ~= 3
        error('Incorrect dimensions of the image')
    end
    if size(img, 1) ~= size(img, 2) || size(img, 2) ~= size(img, 3) || size(img, 1) ~= size(img, 3)
        error('The input image should have same size of each dimension')
    end
end

if size(varargin, 2) == 1
    if size(varargin{1}, 2) == 3
        order = 2;
        poi = varargin{1};
        if any(poi < 0) || any(poi > size(img, 1))
            error('Incorrect value of the point-of-interest')
        end
    else
        poi = [(size(img,1)-1)/2, (size(img,2)-1)/2, (size(img,3)-1)/2];
        order = varargin{1};
        if order < 0
            error('Incorrect value of the order')
        end
    end
end

if size(varargin, 2) == 2
    if size(varargin{1}, 2) == 3
        poi = varargin{1};
        order = varargin{2};
    else
        order = varargin{1};
        poi = varargin{2};
    end
end

const = prepStep(size(img, 1), order);
invariants = extractFeatures(img, poi(1), poi(2), poi(3), const);
    
    
