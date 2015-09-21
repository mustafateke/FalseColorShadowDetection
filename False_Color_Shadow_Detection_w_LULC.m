% LULC
% Shadow Detection from False Color NIR Images
% Author: Mustafa Teke, mustafa.teke@gmail.com
% Ankara, Turkey, 2014

% Notes: 
% Results with raw data (Digital Number) vs Radiometrically corrected 
% Data (Radiance/Reflectance) may vary. 
% For example Landsat data gives better results with radiance/reflectance
% If you are working with Digital Numbers, you may want to stretch data.
% Also note that if you are working with Pansharpened image, some
% pansharpeining algorithms change radiometry of an image. 

% If you need, you may use morphological operations for high resolution
% data, this may not need for Landsat

% If you use this method please cite the following article
% Teke, M.; Baþeski, E.; Ok, A. Ö.; Yüksel, B. & Þenaras, Ç., 
% Multi-spectral false color shadow detection Photogrammetric Image Analysis, 
% Springer Berlin Heidelberg, 2011, 109-119

clear all;
close all;

filename='LC81800312014262LGN00_ref_sub.tif';
im = imread(filename);
imgI=im(:,:,5);
imgR=im(:,:,4);
imgG=im(:,:,3);
imgB=im(:,:,2);


%% Local Stretching
% maxImValue=2048;
% imgR=double(imgR-minR);
% imgR=double(imgR)/max(imgR(:)) ;
% imgR=uint16(imgR*maxImValue);
% 
% imgG=double(imgG-min(imgG(:))) ;
% imgG=imgG/max(imgG(:)) ;
% imgG=uint16(imgG*maxImValue);
% 
% imgB=double(imgB-min(imgB(:)));
% imgB=imgB/max(imgB(:)) ;
% imgB=uint16(imgB*maxImValue);
% 
% imgI=double(imgI-min(imgI(:)));
% imgI=imgI/max(imgI(:)) ;
% imgI=uint16(imgI*maxImValue);

%% Global Stretching
minR = min(imgR(:));
minG = min(imgG(:));
minB = min(imgB(:));
minI = min(imgI(:));
maxR = max(imgR(:));
maxG = max(imgG(:));
maxB = max(imgB(:));
maxI = max(imgI(:));

minVal = double(min([minR minG minB minI]));
maxVal = double(max([maxR maxG maxB maxI]));

imgR=double(imgR-minVal);
imgR=imgR/maxVal ;
% imgR=uint16(imgR*maxImValue);

imgG=double(imgG-minVal) ;
imgG=imgG/maxVal ;
% imgG=uint16(imgG*maxImValue);

imgB=double(imgB-minVal);
imgB=imgB/maxVal ;
% imgB=uint16(imgB*maxImValue);

imgI=double(imgI-minVal);
imgI=imgI/maxVal ;
% imgI=uint16(imgI*maxImValue);

cloudThresh = 0.45;
%% Cloud Map
CloudMapB = im2bw(imadjust(imgB), cloudThresh); %figure, imshow(CloudMapB);
CloudMapG = im2bw(imadjust(imgG), cloudThresh); %figure, imshow(CloudMapG);
CloudMapR = im2bw(imadjust(imgR), cloudThresh); %figure, imshow(CloudMapR);
CloudMapI = im2bw(imadjust(imgI), cloudThresh); %figure, imshow(CloudMapI);

CloudMap = ...
    CloudMapB & ...
    CloudMapG & ...
    CloudMapR & ...
    CloudMapI ...
    ;
figure, imshow(CloudMap);
%% Save RGB Image
imgRGB(:,:,1) = histeq(imgR);
imgRGB(:,:,2) = histeq(imgG);
imgRGB(:,:,3) = histeq(imgB);
imwrite(imgRGB, [filename '_RGB.png']);
figure, imshow(imgRGB);

%% Save False Color Image
imgIRG(:,:,1) = imadjust(imgI);
imgIRG(:,:,2) = imadjust(imgR);
imgIRG(:,:,3) = imadjust(imgG);
imwrite(imgIRG, [filename '_IRG.png']);
figure, imshow(imgIRG);

%% Shadow Detection
% For optimal shadow detection imadjust, saturating %1 data works optimal
% However imadjusting data for shadow and water does not produce optimal
% results
r=imadjust((double(imgI)));
g=imadjust((double(imgR)));
b=imadjust((double(imgG)));
[h,s,v] = rgb2hsv(r,g,b);
ShadowIndex = (s - v)./(v+s);
figure; imshow(ShadowIndex) ; title('Raw Shadow Map');
imwrite(uint8(127*(ShadowIndex+1)), [filename '_Shadow_Index.png']);

%% Water Detection
NDWI = (double(imgG)-double(imgI))./(double(imgG)+double(imgI));
imwrite(uint8(127*(NDWI+1)), [filename '_NDWI.png']);
NDWI(imgI == 0) = 0;
figure; imshow(NDWI) ; title('Raw Water Map');

waterthresh = graythresh(NDWI);
WaterMap = im2bw(NDWI, waterthresh);
figure, imshow(WaterMap);title('Thresholded Water Map');

%% Vegetation Detection
NDVI=(double(imgI)-double(imgR))./(double(imgI)+double(imgR));
imwrite(uint8(127*(NDVI+1)), [filename '_NDVI.png']);
NDVI(imgI == 0) = 0;
figure; imshow(NDVI) ; title('Raw Vegetation Map');

vegetationthresh = graythresh(NDVI);
VegetationMap = im2bw(NDVI, vegetationthresh);
vegetationthresh2 = graythresh(NDVI(VegetationMap>0));
% this second thresholding is needed to differentiate shadowed vegetation areas.
% Comment if your data is different
VegetationMap2 = im2bw(NDVI, vegetationthresh2); 

figure, imshow(VegetationMap);title('Thresholded Vegetation Map');
figure, imshow(VegetationMap2);title('Thresholded Vegetation Map2');

%% Final Shadow Map
shadowthresh = graythresh(ShadowIndex);
ShadowMap = im2bw(ShadowIndex, shadowthresh);
figure; imshow(ShadowMap) ; title('Thresholded Shadow Map');

FinalShadowMap = 255*uint8(ShadowMap-VegetationMap2-WaterMap);
figure, imshow(uint8(FinalShadowMap));title('Final Shadow Map');  
imwrite(uint8(FinalShadowMap), [filename '_Shadow_Map.png']);

[rows, cols, bands] = size(FinalShadowMap);

LULCMap_Red = 165*ones(rows, cols);
LULCMap_Green = 42*ones(rows, cols);
LULCMap_Blue = 42*ones(rows, cols);

LULCMap_Red(VegetationMap > 0) = 0;
LULCMap_Green(VegetationMap > 0) = 255;
LULCMap_Blue(VegetationMap > 0) = 0;

LULCMap_Red(WaterMap > 0) = 0;
LULCMap_Green(WaterMap > 0) = 0;
LULCMap_Blue(WaterMap > 0) = 255;

LULCMap_Red(CloudMap > 0) = 255;
LULCMap_Green(CloudMap > 0) = 255;
LULCMap_Blue(CloudMap > 0) = 255;

LULCMap_Red(FinalShadowMap > 0) = 0;
LULCMap_Green(FinalShadowMap > 0) = 0;
LULCMap_Blue(FinalShadowMap > 0) = 0;

LULCMap(:,:,1) = LULCMap_Red;
LULCMap(:,:,2) = LULCMap_Green;
LULCMap(:,:,3) = LULCMap_Blue;
figure, imshow(uint8(LULCMap));
imwrite(uint8(LULCMap), [filename '_LULC.png']);

ShadowPercent = 100*length(FinalShadowMap(FinalShadowMap>0))/(rows*cols)
CloudPercent = 100*length(CloudMap(CloudMap>0))/(rows*cols)
