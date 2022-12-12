clear all; close all;

variables = load("connectionOutput.mat");
originalImage = variables.originalImage;
boxImage = variables.boxImage;
variances = variables.variances;
meanVariances = variables.meanVariances;
localWindowVariancesMatrix = variables.localConnectedVariancesMatrix;

figure; imshow(variables.boxImage, "InitialMagnification", 'fit'); title("Noisy Image");

image = boxImage;
gaussianImage = image;


disp("Size="+size(image));
meanValue = mean(image, 'all');
varianceValue = var(double(image), 1, 'all');
disp("Mean="+meanValue);
disp("Variance="+varianceValue);

% adaptive local filtering

disp("Image Size "+size(boxImage,1)+" "+size(boxImage,2));
image = double(gaussianImage);

noislessImage1 = adaptiveLocalNoiseReductionWithNoiseVariance(image);
noislessImage2 = adaptiveLocalNoiseReductionWithGlobalNoiseVariance(image);
noislessImage3 = adaptiveLocalNoiseReductionWithMeanLocalWindowNoiseVariance(image);
noislessImage4 = adaptiveLocalNoiseReductionWithConnectedNoiseVariance(image, variances, localWindowVariancesMatrix);
noislessImage5 = adaptiveLocalNoiseReductionWithConnectedNoiseVariance2(image, variances, meanVariances, localWindowVariancesMatrix);

figure; imshow(originalImage, "InitialMagnification", 'fit'); title("Original Image");
figure; imshow(boxImage, "InitialMagnification", 'fit'); title("Noisy Image");


figure; imshow(uint8(noislessImage1), "InitialMagnification", 'fit'); title("ALF - Noise variance");
figure; imhist(noislessImage1, 255);title("ALF - Noise variance PDF");

figure; imshow(uint8(noislessImage2), "InitialMagnification", 'fit'); title("ALF - Global variance");
figure; imhist(noislessImage2, 255);title("ALF - Global variance PDF");

figure; imshow(uint8(noislessImage3), "InitialMagnification", 'fit'); title("ALF - Mean local window variances");
figure; imhist(noislessImage3, 255);title("ALF - Mean local window variances PDF");

figure; imshow(uint8(noislessImage4), "InitialMagnification", 'fit'); title("ALF - Mean connected local window variances");
figure; imhist(noislessImage4, 255);title("ALF - Mean local window variances PDF");

figure; imshow(uint8(noislessImage5), "InitialMagnification", 'fit'); title("ALF - Connected adaptable local window variances");
figure; imhist(noislessImage5, 255);title("ALF - Connected adaptable local window variances PDF");

save("adaptiveOutput");

function noislessImage1 = adaptiveLocalNoiseReductionWithNoiseVariance(image)
    % assuming we know the noise variance
    noiseVariance = 0.001;
    
    % defining window
    windowX = 3;
    windowY = 3;
    
    paddedImage = padarray(image, [1 1], 0);

    noislessImage = zeros(size(image));
    
    [localMeanMatrix, localVarianceMatrix] = imageStatistics(image);
    
    % running the window over the image
    for i = 1:(size(paddedImage,1)-windowX+1)
        for j = 1:(size(paddedImage,2)-windowY+1)
            windowSelection = paddedImage(i:i+windowX-1, j:j+windowY-1);
    
            localMeanMatrix(i,j) = mean(windowSelection, "all");
            localVarianceMatrix(i,j) = var(windowSelection, 1, "all");
            disp("LV="+localVarianceMatrix(i,j)+" NV="+noiseVariance);

            % Noise Variance should be less that or equal to neighbourhood
            noiseVarianceRatio = 1;
            if noiseVariance <= localVarianceMatrix(i,j)
                noiseVarianceRatio = noiseVariance/localVarianceMatrix(i,j);
            end

    
            noislessImage(i,j) = image(i,j) - ((noiseVarianceRatio) * (image(i,j) - localMeanMatrix(i,j)));
    
        end
    end

    noislessImage1 =  uint8(noislessImage);

end

function noislessImage2 = adaptiveLocalNoiseReductionWithGlobalNoiseVariance(image)
    
    % assuming we know the noise variance
    % noise variance is global noise variance of image
    noiseVariance = var(image, 1, "all");
    
    % defining window
    windowX = 3;
    windowY = 3;
    
    paddedImage = padarray(image, [1 1], 0);

    noislessImage = zeros(size(image));

    [localMeanMatrix, localVarianceMatrix] = imageStatistics(image);

    % running the window over the image
    for i = 1:(size(paddedImage,1)-windowX+1)
        for j = 1:(size(paddedImage,2)-windowY+1)
            windowSelection = paddedImage(i:i+windowX-1, j:j+windowY-1);
    
            localMeanMatrix(i,j) = mean(windowSelection, "all");
            localVarianceMatrix(i,j) = var(windowSelection, 1, "all");
            disp("LV="+localVarianceMatrix(i,j)+" NV="+noiseVariance);

            noiseVarianceRatio = 1;
            if noiseVariance <= localVarianceMatrix(i,j)
                noiseVarianceRatio = noiseVariance/localVarianceMatrix(i,j);
            end
    
            noislessImage(i,j) = image(i,j) - ((noiseVarianceRatio) * (image(i,j) - localMeanMatrix(i,j)));
    
        end
    end

    noislessImage2 = uint8(noislessImage);

end

function noislessImage3 = adaptiveLocalNoiseReductionWithMeanLocalWindowNoiseVariance(image)
    
    % defining window
    windowX = 3;
    windowY = 3;
    
    paddedImage = padarray(image, [1 1], 0);

    noislessImage = zeros(size(image));

    [localMeanMatrix, localVarianceMatrix] = imageStatistics(image);
    
    % assuming we know the noise variance
    % noise variance is mean of local window noise variance of image
    noiseVariance = mean(localVarianceMatrix,"all");

    % running the window over the image
    for i = 1:(size(paddedImage,1)-windowX+1)
        for j = 1:(size(paddedImage,2)-windowY+1)
            windowSelection = paddedImage(i:i+windowX-1, j:j+windowY-1);
    
            localMeanMatrix(i,j) = mean(windowSelection, "all");
            localVarianceMatrix(i,j) = var(windowSelection, 1, "all");
            disp("LV="+localVarianceMatrix(i,j)+" NV="+noiseVariance);

            noiseVarianceRatio = 1;
            if noiseVariance <= localVarianceMatrix(i,j)
                noiseVarianceRatio = noiseVariance/localVarianceMatrix(i,j);
            end
    
            noislessImage(i,j) = image(i,j) - ((noiseVarianceRatio) * (image(i,j) - localMeanMatrix(i,j)));
    
        end
    end

    noislessImage3 = uint8(noislessImage);
end

function noislessImage4 = adaptiveLocalNoiseReductionWithConnectedNoiseVariance(image, connectedVariances, localWindowVariancesMatrix)
    % noise variance through connections
    
    
    % defining window
    windowX = 3;
    windowY = 3;
    
    paddedImage = padarray(image, [1 1], 0);

    noislessImage = zeros(size(image));

    [localMeanMatrix, localVarianceMatrix] = imageStatistics(image);

    % running the window over the image
    for i = 1:(size(paddedImage,1)-windowX+1)
        for j = 1:(size(paddedImage,2)-windowY+1)
            % noise variance is variance value of independent connections
            noiseVariance = connectedVariances(i,j);

            windowSelection = paddedImage(i:i+windowX-1, j:j+windowY-1);
    
            localMeanMatrix(i,j) = mean(windowSelection, "all");
            localVarianceMatrix(i,j) = var(windowSelection, 1, "all");
            disp("LV="+localVarianceMatrix(i,j)+" NV="+noiseVariance);

            noiseVarianceRatio = 1;
            if noiseVariance <= localVarianceMatrix(i,j)
                noiseVarianceRatio = noiseVariance/localVarianceMatrix(i,j);
            end
    
            noislessImage(i,j) = image(i,j) - ((noiseVarianceRatio) * (image(i,j) - localMeanMatrix(i,j)));

%             if noiseVariance <= localWindowVariancesMatrix(i,j)
%                 noiseVarianceRatio = noiseVariance/localWindowVariancesMatrix(i,j);
%             end
%     
%             noislessImage(i,j) = image(i,j) - ((noiseVarianceRatio) * (image(i,j) - localMeanMatrix(i,j)));
    
        end
    end

    noislessImage4 = uint8(noislessImage);
end

function noislessImage5 = adaptiveLocalNoiseReductionWithConnectedNoiseVariance2(image, connectedVariances, meanVariances, localWindowVariancesMatrix)
    % noise variance through connections
    
    
    % defining window
    windowX = 3;
    windowY = 3;
    
    paddedImage = padarray(image, [1 1], 0);

    noislessImage = zeros(size(image));

    [localMeanMatrix, localVarianceMatrix] = imageStatistics(image);

    % running the window over the image
    for i = 1:(size(paddedImage,1)-windowX+1)
        for j = 1:(size(paddedImage,2)-windowY+1)
            % noise variance is variance value of independent connections
            noiseVariance = connectedVariances(i,j);
            noiseVariance = (meanVariances(i,j)+connectedVariances(i,j));

            windowSelection = paddedImage(i:i+windowX-1, j:j+windowY-1);
    
            localMeanMatrix(i,j) = mean(windowSelection, "all");
            localVarianceMatrix(i,j) = var(windowSelection, 1, "all");
            disp("LV="+localVarianceMatrix(i,j)+" NV="+noiseVariance);

            noiseVarianceRatio = 1;

            localVariance = (localVarianceMatrix(i,j));
            if noiseVariance <= localVariance
                noiseVarianceRatio = noiseVariance/localWindowVariancesMatrix(i,j);
                noiseVarianceRatio = noiseVariance/localVariance;
            end
    
            noislessImage(i,j) = image(i,j) - ((noiseVarianceRatio) * (image(i,j) - localMeanMatrix(i,j)));
    
        end
    end

    noislessImage5 = uint8(noislessImage);
end

function [localMeanMatrix, localVarianceMatrix] = imageStatistics(image)


    localMeanMatrix = zeros(size(image));
    localVarianceMatrix = zeros(size(image));

    windowX = 3;
    windowY = 3;
    
    paddedImage = padarray(image, [1 1], 0);

    % running the window over the image
    for i = 1:(size(paddedImage,1)-windowX+1)
        for j = 1:(size(paddedImage,2)-windowY+1)
            windowSelection = paddedImage(i:i+windowX-1, j:j+windowY-1);
    
            localMeanMatrix(i,j) = mean(windowSelection, "all");
            localVarianceMatrix(i,j) = var(windowSelection, 1, "all");
            
        end
    end
end






