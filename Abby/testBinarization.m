pathToWatch = 'C:\Users\Abby\Cornell\Abby_test_tethered_code\01_23082024\' ; 
allBGcell = findBGTethered(pathToWatch);

LoadPhantomLibraries();
RegisterPhantom(true);

movNum = 1;
movNumStr = num2str(movNum,'%03.f');
cineSuffix = ['_',movNumStr,'.cine'];
movieNum = movNumStr;

camNames = {'yz','xz','xy'};
camFilenames = strcat(camNames,cineSuffix);
cinFilenames = fullfile(pathToWatch,camFilenames);

camInd = 3;

cineData = myOpenCinFile(cinFilenames{camInd});
currCineIm = myReadCinImage(cineData,0);
bg = allBGcell{camInd};

subtractIm = bg-currCineIm;

% x_COM = [300,320,320];
% y_COM = [340,350,400];
% windowLength = 60;
% 
% mask = false(size(bg));
% mask(y_COM(camInd)-windowLength:y_COM(camInd)+windowLength,...
%     x_COM(camInd)-windowLength:x_COM(camInd)+windowLength) = true;
% 
% test1 = subtractIm;
% test1(~mask) = false;

mask=imbinarize(subtractIm);
imshow(mask)

% imshow(test1)
% level = graythresh(test1)*0.6;
% imshow(subtractIm)
% figure()
% imshow(imbinarize(test1,level))
UnregisterPhantom();
UnloadPhantomLibraries();


