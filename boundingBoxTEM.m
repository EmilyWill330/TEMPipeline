function extractedNP = boundingBoxTEM(properties, IAdj)
% Function used in the final pipeline to output TEM images of individual
% nanoparticles for evaluation by a trained neural network
% Coded by Aaron Ghrist
% Updated by Aaron Ghrist on 8-29-2021

extractedNP = uint8.empty(44,44,0,length(properties.Area));

for i = 1:length(properties.Area) %Cycles through identified (non-overlapping) particles
    
    % Need the bounding boxes for each to export image
    index1 = properties.BoundingBox(i,1:2) - [2,2];
    index2 = index1 + properties.BoundingBox(i,3:4) + [4,4];
    % Pads indices of bounding boxes by two pixels in every direction to avoid cropping boundary information
    
    % Boundary Checks
    if floor(index1(1)) < 1
        x1 = 1;
    else
        x1 = floor(index1(1));
    end
    if floor(index1(2)) < 1
        y1 = 1;
    else
        y1 = floor(index1(2));
    end
    if ceil(index2(1)) > size(IAdj,1)
        x2 = size(IAdj,2);
    else
        x2 = ceil(index2(1));
    end
    if ceil(index2(2)) > size(IAdj,2)
        y2 = size(IAdj,1);
    else
        y2 = ceil(index2(2));
    end
    
    % Extract relevant region of pre-processed TEM image
    cropped = IAdj(x1:x2, y1:y2);
    
    % Obtains current aspect ratio of image (height divided by width)
    aspectRatio = length(cropped(:,1))/ length(cropped(1,:));
    
    % Standardizes by downscaling to match simulated images
    imageScaled = uint8.empty(44,44,0); %length(properties.Area)
    
    % Maintains the proper aspect ratio while converting the entire image to a square
    if aspectRatio >= 1
        % Pads the image with zeros to maintain the ratio
        resizeTemp = imresize(cropped, [44, round(44/aspectRatio)]);
        imageScaled(:,(23-floor(length(resizeTemp(1,:))/2)):(22-floor(length(resizeTemp(1,:))/2))+length(resizeTemp(1,:)),1) = resizeTemp;
    else
        resizeTemp = imresize(cropped, [round(44*aspectRatio),44]);
        resizeTemp = rot90(resizeTemp);
        imageScaled(:,(23-floor(length(resizeTemp(1,:))/2)):(22-floor(length(resizeTemp(1,:))/2))+length(resizeTemp(1,:)),1) = resizeTemp;
    end
    
    % Stochastically corrects for the black bars present when aspectRatio is not 1 (helps significantly with normalization later)
    if aspectRatio ~= 1
        leftpts = imageScaled(:,(23-floor(length(resizeTemp(1,:))/2)):(27-floor(length(resizeTemp(1,:))/2)));
        imageScaled(:,1:(22-floor(length(resizeTemp(1,:))/2))) = round(0.5*mean(leftpts(:))) + 10*poissrnd(5, 44, (22-floor(length(resizeTemp(1,:))/2)));
        rightpts = imageScaled(:,(18-floor(length(resizeTemp(1,:))/2))+length(resizeTemp(1,:)):(22-floor(length(resizeTemp(1,:))/2))+length(resizeTemp(1,:)));
        imageScaled(:,(23-floor(length(resizeTemp(1,:))/2))+length(resizeTemp(1,:)):end) = round(0.5*mean(rightpts(:))) + 10*poissrnd(5, 44, (45-((23-floor(length(resizeTemp(1,:))/2))+length(resizeTemp(1,:)))));
    end
    
    % Randomly rotates the matrix 90 degrees 1-4 times to correct for the possible systemic error introduced with aspect-ratio padding
    imageScaled = rot90(imageScaled, randi(4));
    
    % Normalizes the data for optimal NN evaluation (no need to invert color scheme like in simulation; this is already correct)
    maxIntensity = max(imageScaled(:)); minIntensity = min(imageScaled(:));
    
    % Exports to the 4D array
    extractedNP(:,:,1,i) = (imageScaled-minIntensity)*(255/(maxIntensity-minIntensity)); 
end

end