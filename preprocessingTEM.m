function [properties, IAdj, particleSep, particleOverlay] = preprocessingTEM(image)
% Function used in the final pipeline to quickly process TEM images for
% extracting properties
% Coded by Emily Williamson
% Updated by Aaron Ghrist on 8-29-2021

I = im2gray(image);
IAdj = imadjust(I);

% Filtering
F = adapthisteq(IAdj,"Range","original");

F = medfilt2(F,[10 10],"symmetric");

% show = labeloverlay(Ifilt,seg); % Change this threshold depending on material
IAdjBW = imbinarize(F);

SE1 = strel("disk",5);
IAdjBWc = imdilate(IAdjBW,SE1);

IAdjopen = bwareaopen(~IAdjBWc,900);

IAdjC = imclearborder(IAdjopen,4);

BW_fill_filter = imdilate(IAdjC,SE1);

% imFill = imfill(imDil,"holes");
% 
% BW_fill_filter = imclose(imFill,SE1);

%regions = bwconncomp(BW_fill_filter);
%labels = labelmatrix(regions);

%labelsRGB = label2rgb(labels);

dist = bwdist(~BW_fill_filter);
dist = imcomplement(dist);

dist = imhmin(dist,1);
particleSep = watershed(dist);
particleSep(~BW_fill_filter) = 0;

particleOverlay = labeloverlay(IAdj,particleSep);

properties = regionprops("table",particleSep,"all");
% keep = properties_full.Solidity > 0.90;
% properties = properties_full(keep,:);


end