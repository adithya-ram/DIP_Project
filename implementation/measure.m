clear all; close all;

connectionVariables = load("connectionOutput.mat");
connectivityMatrices = connectionVariables.connectionsMatrix;
connectionsCount = connectionVariables.connectionsCount;
adaptiveVariables = load("adaptiveOutput.mat");

originalImage = connectionVariables.originalImage;
boxImage = connectionVariables.boxImage;

noiselessImage1 = adaptiveVariables.noislessImage1;
noiselessImage2 = adaptiveVariables.noislessImage2;
noiselessImage3 = adaptiveVariables.noislessImage3;
noiselessImage4 = adaptiveVariables.noislessImage4;
noiselessImage5 = adaptiveVariables.noislessImage5;

% Mean Filter
w = fspecial('average',3);
averagedNoiseImage= conv2(boxImage, w, 'same');
averagedNoiseImage=uint8(averagedNoiseImage);

% Median Filter
medianedNoiseImage = medfilt2(boxImage,[3,3]);
medianedNoiseImage=uint8(medianedNoiseImage);

diff = abs(originalImage - boxImage);

diffMean = abs(originalImage - averagedNoiseImage);
diffMedian = abs(originalImage - medianedNoiseImage);

diff1 = abs(originalImage - noiselessImage1);
diff2 = abs(originalImage - noiselessImage2);
diff3 = abs(originalImage - noiselessImage3);
diff4 = abs(originalImage - noiselessImage4);
diff5 = abs(originalImage - noiselessImage5);

totalNoiseError = sum(diff, 'all');
error1 = sum(diff1, 'all');
error2 = sum(diff2, 'all');
error3 = sum(diff3, 'all');
error4 = sum(diff4, 'all');
error5 = sum(diff5, 'all');

errorMean = sum(diffMean, 'all');
errorMedian = sum(diffMedian, 'all');

disp("Noisless Main Errors: "+totalNoiseError);
disp("Noisless Image 1 Errors: "+error1);
disp("Noisless Image 2 Errors: "+error2);
disp("Noisless Image 3 Errors: "+error3);
disp("Noisless Image 4 Errors: "+error4);
disp("Noisless Image 5 Errors: "+error5);
disp("Noisless Mean Errors: "+errorMean);
disp("Noisless Median Errors: "+errorMedian);

per1 = (error1/totalNoiseError) * 100;
per2 = (error2/totalNoiseError) * 100;
per3 = (error3/totalNoiseError) * 100;
per4 = (error4/totalNoiseError) * 100;
per5 = (error5/totalNoiseError) * 100;

perMean = (errorMean/totalNoiseError) * 100;
perMedian = (errorMedian/totalNoiseError) * 100;

disp("Noisless Image 1 percentage: "+per1);
disp("Noisless Image 2 percentage: "+per2);
disp("Noisless Image 3 percentage: "+per3);
disp("Noisless Image 4 percentage: "+per4);
disp("Noisless Image 5 percentage: "+per5);

disp("Noisless Mean percentage: "+perMean);
disp("Noisless Median percentage: "+perMedian);

