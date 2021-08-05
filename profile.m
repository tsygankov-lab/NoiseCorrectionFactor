%% Dr. Denis Tsygankov (2021)
%% Integrative Systems Biology Lab, Wallace H. Coulter Department of Biomedical Engineering, Georgia Institute of Technology and Emory University SOM 

%% If you use any part of this script, please cite:
%
% DJ Marston et al. "Correcting artifacts in ratiometric biosensor imaging; an improved approach for dividing noisy signals"
% Frontiers in Cell and Developmental Biology (2021), doi: 10.3389/fcell.2021.685825

%% This function calculates mean ratio signal along level-set contours as a function of the distance from the cell edge 

%% Required Arguments:
%
%  im - intensity image to process
%  mask - cell mask image (must have the same dimensions as im)
%  Depth - the number of pixels to go from the cell edge
%  inout - specifies if the calculation is for inside or outside the cell: 
%          1 sets the direction from the cell edge toward the cell center
%          0 sets the direction from the cell edge outward the cell center

%% Output: (see examples of the output in Figures 1D,1E,3,4E,6C,6D of the paper)
%
%  x - x-coordinates of the calculated mean ratio function 
%  means - y-coordinates of the calculated mean ratio function
%  meansAddStd - mean plus standard deviation
%  meansSubStd  - mean minus standard deviation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,means,meansAddStd,meansSubStd] = profile(im,mask,Depth,inout) 
    
    meansVal = cell(1,Depth);
    stdsVal = cell(1,Depth);
    interShrink = mask;
    interShrink = imfill(interShrink, 'holes');
    
    if inout
        boundaryMask = bwperim(interShrink);
    else        
        boundaryMask = bwperim(~interShrink);
        boundaryMask(1,1:end) = 0;
        boundaryMask(end,1:end) = 0;
        boundaryMask(1:end,1) = 0;
        boundaryMask(1:end,end) = 0;
    end
    
    for i = 1:Depth
        pixels = double(im(boundaryMask));
        meansVal{i} = mean(pixels(pixels~=0));
        stdsVal{i} = std(pixels(pixels~=0));
        
        if inout
            interShrink = imerode(mask, strel('disk', i));
            boundaryMask = bwperim(interShrink);
        else
            interShrink = imdilate(mask, strel('disk', i));
            boundaryMask = bwperim(~interShrink);
            boundaryMask(1,1:end) = 0;
            boundaryMask(end,1:end) = 0;
            boundaryMask(1:end,1) = 0;
            boundaryMask(1:end,end) = 0;
        end
    end
    means = [meansVal{:}];
    stds = [stdsVal{:}];
    meansAddStd = means + stds;
    meansSubStd = means - stds;
    
    if inout
        x = 0:(Depth-1);
    else
        x = -1:(-1):(-Depth);
    end
  
end