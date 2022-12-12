clear all; close all;

boxImage = zeros(100,100);
boxImage = uint8(boxImage+100);
boxImage(25:75, 25:75) = 150;

boxImage = zeros(100,100);
boxImage = uint8(boxImage+100);
boxImage(20:80, 20:80) = 150;
boxImage(40:60, 40:60) = 100;


% boxImage = zeros(10,10);
% boxImage = uint8(boxImage+100);
% boxImage(4:7, 4:7) = 150;

% boxImage = zeros(7,7);
% boxImage = uint8(boxImage+100);
% boxImage(3:5, 3:5) = 150;

originalImage = boxImage;
figure; imshow(originalImage, "InitialMagnification", 'fit'); title("Original Image");

figure; imhist(originalImage, 255);title("Original Image PDF");

% boxImage = imread("../ImageBank/testpattern512.tif");
% boxImage = imread("../ImageBank/lena.tif");

% noising the image
boxImage = imnoise(boxImage, "gaussian", 0, 0.001);
% boxImage = imnoise(boxImage, "salt & pepper", 0.01);
image = boxImage;

% boxImage = zeros(10,10);
% boxImage = uint8(boxImage+100);
% boxImage(4:7, 4:7) = 150;
% % noising the image
% boxImage = imnoise(boxImage, "gaussian", 0, 0.001);
% % boxImage = imnoise(boxImage, "salt & pepper", 0.1);

% noised Image
figure;imshow(boxImage, "InitialMagnification", 'fit'); title("Noised Image"); 
figure; imhist(boxImage, 255);title("Noised Image PDF");

% after noising, run mean or median filterning to get the connections
% better
% Kernel=3
% w= fspecial('average',Kernel) 
% averagedImage= conv2(boxImage, w, 'same');
% averagedImage=uint8(averagedImage);
% image = averagedImage;
% figure;imshow(image, "InitialMagnification", 'fit'); title("Unnoised Image");

medianedImage = medfilt2(boxImage,[3,3]);
medianedImage=uint8(medianedImage);
image = medianedImage;
figure;imshow(image, "InitialMagnification", 'fit'); title("Unnoised Image for 4-connectivity");

% boxImage = zeros(4,4);
% boxImage = uint8(boxImage+100);
% boxImage(2:3, 2:3) = 150;
% % noising the image
% boxImage = imnoise(boxImage, "gaussian", 0, 0.001);
% image = boxImage;



figure; imshow(image, "InitialMagnification", 'fit'); title("Unnoised image");

% figure; imhist(image);
image = image > 125;
figure; imshow(image, "InitialMagnification", 'fit'); title("Logical Unnoised image");

[m,n] = size(image);

flag = true;


global travelledMatrix;
travelledMatrix = zeros(size(boxImage));
travelledMatrixCoverage = find(travelledMatrix==0);

global connectionsCount;
connectionsCount = 1;
global conMatrix;
conMatrix = zeros([m,n,10]);
global connectionsMatrix;
connectionsMatrix = zeros([m,n,10]);


while(~isempty(travelledMatrixCoverage))
    [row, column, data] = find(travelledMatrix == 0);
    connectionStartX = row(1);
    connectionStartY = column(1);
    
    conMatrix(connectionStartX,connectionStartY,connectionsCount) = 1;

    tic
    connectionsMatrix(:,:,connectionsCount) = connectivity(image, conMatrix(:,:,connectionsCount), travelledMatrix, connectionStartX, connectionStartY);
    timeElapsed = toc
    
    connectionsCount = connectionsCount + 1;
    travelledMatrixCoverage = find(travelledMatrix==0);
end


figure; 
subplot(1,2,1); imshow(boxImage); title("Original Image");
subplot(1,2,2); imshow(image); title("Logical Image");

displayConnections(connectionsMatrix, connectionsCount-1, boxImage);

variances = captureConnectionsVariances(connectionsMatrix, connectionsCount-1, boxImage);
[localConnectedVariancesMatrix, meanVariances] = captureConnectionsLocalVariances(connectionsMatrix, connectionsCount-1, boxImage);

save("connectionOutput")

function [connectionsMatrix] = connectivity(image, cm, tm, i, j)
    [m,n] = size(image);
    connectionsMatrix = cm
    if i<1 || i>m || i<1 || i>m || j <1 || j>n || j<1 || j>n
        return
    end

    global travelledMatrix;
    travelledMatrix = tm;
    
    if travelledMatrix(i,j) == 1
        return
    else
        travelledMatrix(i,j) = 1;
    end

    disp("Image iteration i="+i+" j="+j);
    currentIndexValue = image(i,j);

    
    
    top = i-1;
    bottom = i+1;
    left = j-1;
    right = j+1;

    if top>=1 && top <=m
        topIndexValue = image(top,j);
        

        if currentIndexValue == topIndexValue
            cm(top,j) = 1;
            
            if travelledMatrix(top,j) ~= 1
                cm = connectivity(image,cm,travelledMatrix, top,j);
            end
        end
    end

    if bottom >=1 && bottom<=m
        bottomIndexValue = image(bottom,j);
        

        if currentIndexValue == bottomIndexValue
            cm(bottom,j) = 1;
            
            if travelledMatrix(bottom,j) ~= 1
                cm = connectivity(image,cm,travelledMatrix, bottom,j);
            end
        end
    end

    if left>=1 && left <=n
        leftIndexValue = image(i,left);
        

        if currentIndexValue == leftIndexValue
            cm(i,left) = 1;

            if travelledMatrix(i,left) ~= 1
                cm = connectivity(image,cm,travelledMatrix, i,left);
            end
        end
    end

    if right >=1 && right<=n
        rightIndexValue = image(i,right);
        

        if currentIndexValue == rightIndexValue
            cm(i,right) = 1;

            if travelledMatrix(i,right) ~= 1
                cm = connectivity(image,cm,travelledMatrix, i,right);
            end
        end
    end

    connectionsMatrix = cm;
end

function display = displayConnections(connectivityMatrices, connectionsCount, image)

    [m,n] = size(image);
    for connectionIndex = 1:connectionsCount
        connectionMatrix = connectivityMatrices(:,:,connectionIndex);

        mask = connectionMatrix;
        overlapIndices = find(mask==1);

        reshapedImage = reshape(image, [m*n,1]);

        outputImage = zeros(size(image));
        reshapedOutputImage = reshape(outputImage, [m*n,1]);
        reshapedOutputImage(overlapIndices,:) = reshapedImage(overlapIndices,:);
        outputImage = uint8(reshape(reshapedOutputImage, [m,n]));

        figure; 
        imshow(outputImage, "InitialMagnification", 'fit'); title("Connection "+connectionIndex);

    end
    
    display = 1
end

function variances = captureConnectionsVariances(connectivityMatrices, connectionsCount, image)

    [m,n] = size(image);
    variances = zeros(size(image));
    variances(1:m, 1:n) = -1;
    
    
    for connectionIndex = 1:connectionsCount
        connectionMatrix = connectivityMatrices(:,:,connectionIndex);

        mask = connectionMatrix;
        overlapIndices = find(mask==1);

        connectedImage = uint8(mask).*image;
        reshapedConnectedImage = reshape(connectedImage, [m*n, 1]);
        reshapedConnectionPixelValues = reshapedConnectedImage(overlapIndices,:);
        varianceConnectionPixelValues = var(double(reshapedConnectionPixelValues), 1, 'all');

        reshapedVarianceMatrix = reshape(variances, [m*n,1]);
        reshapedVarianceMatrix(overlapIndices,:) = varianceConnectionPixelValues;
        variances = reshape(reshapedVarianceMatrix, [m,n]);

    end
end

% plan is to calculate variances also as a group of 3
% earlier it was whole connected variance, not window
% now its to take the connectedness with window 3
% after each connection iteration, within the inner 2 for loops
% maintain local variable for connection specific variances matrix
% then calculat mean, assign to position of the connectivity mask
% this will be new variances, signal noise
% the local neighborhood variance would be the normal, the unmeaned one
function [localConnectedVariancesMatrix, meanConnectedVariances] = captureConnectionsLocalVariances(connectivityMatrices, connectionsCount, image)
    [m,n] = size(image);
    localConnectedVariancesMatrix = zeros(size(image));
    localConnectedVariancesMatrix(1:m, 1:n) = -1;

    localVariancesMatrix = zeros(size(image));
    localVariancesMatrix(1:m, 1:n) = -1;

    meanConnectedVariances = zeros(size(image));
    meanConnectedVariances(1:m, 1:n) = -1;

    paddedImage = padarray(image, [1 1], 0);
    
    % defining window
    windowX = 3;
    windowY = 3;

    for i = 1:(size(paddedImage,1)-windowX+1)
        for j = 1:(size(paddedImage,2)-windowY+1)

            % windowSelection = paddedImage(i:i+windowX-1, j:j+windowY-1);
            
            windowSelection = [
                paddedImage(i,j), paddedImage(i, j+1), paddedImage(i,j+2);
                paddedImage(i+1,j), paddedImage(i+1, j+1), paddedImage(i+1,j+2);
                paddedImage(i+2,j), paddedImage(i+2, j+1), paddedImage(i+2,j+2);
            ];

            windowSelectionReshaped = reshape(windowSelection, [size(windowSelection,1)*size(windowSelection,2), 1]);

            localVariancesMatrix(i,j) = var(double(windowSelectionReshaped), 1, 'all')'
    
        end
    end
    
    for connectionIndex = 1:connectionsCount
        connectionMatrix = connectivityMatrices(:,:,connectionIndex);

        mask = connectionMatrix;
        overlapIndices = find(mask==1);
        
        connectedImage = uint8(mask).*image;
        paddedImage = padarray(connectedImage, [1 1], 0);
    
        % running the window over the image
        for i = 1:(size(paddedImage,1)-windowX+1)
            for j = 1:(size(paddedImage,2)-windowY+1)

                if(mask(i,j)==1)

                    % windowSelection = paddedImage(i:i+windowX-1, j:j+windowY-1);
                
                    windowSelection = [
                        paddedImage(i,j), paddedImage(i, j+1), paddedImage(i,j+2);
                        paddedImage(i+1,j), paddedImage(i+1, j+1), paddedImage(i+1,j+2);
                        paddedImage(i+2,j), paddedImage(i+2, j+1), paddedImage(i+2,j+2);
                    ];
    
                    windowSelectionReshaped = reshape(windowSelection, [size(windowSelection,1)*size(windowSelection,2), 1]);
                    
    
                    windowMask = find(windowSelection ~= 0);
                    windowConnectedPixelValues = windowSelectionReshaped(windowMask,:);
                    windowConnectedVarianceValues = var(double(windowConnectedPixelValues), 1, 'all');
    
                    if ~isnan(windowConnectedVarianceValues)
                        if windowConnectedVarianceValues ~= 0
                            localConnectedVariancesMatrix(i,j) = windowConnectedVarianceValues;
                            meanConnectedVariances(i,j) = windowConnectedVarianceValues;
                        end
                    end

                end
            end
        end

        reshapedMeanConnectedVariances = reshape(meanConnectedVariances, [m*n, 1]);
        varianceMean = mean(reshapedMeanConnectedVariances(overlapIndices,:), 'all');
        reshapedMeanConnectedVariances(overlapIndices,:) = varianceMean;
        meanConnectedVariances = reshape(reshapedMeanConnectedVariances, [m,n]);

    end
end
