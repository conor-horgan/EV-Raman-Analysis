% -------------------------------------------------------------------------
%  Name: EV_Image_Analysis.m
%  Version: 1.0
%  Environment: Matlab 2016a
%  Date: 19/07/2019
%  Author: Conor Horgan
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%  1. Select files for analysis
% -------------------------------------------------------------------------
% Select whole cell image to create mask
uiwait(msgbox('Select a whole cell image to import for image analysis and processing', 'Select Whole Cell Image', 'modal'));
[Cell_Image_File, Cell_Image_Path] = uigetfile({'*.fig;*.png;*.jpg;*.bmp',...
    'Figures (*.fig;*.png;*.jpg;*.bmp)'; '*.*','All Files (*.*)'},'Select Whole Cell Image');

% Select corresponding deuterium image
uiwait(msgbox('Select the corresponding deuterium image to import for image analysis and processing', 'Select Deuterium Image', 'modal'));
[Deuterium_Image_File, Deuterium_Image_Path] = uigetfile({'*.fig;*.png;*.jpg;*.bmp',...
    'Figures (*.fig;*.png;*.jpg;*.bmp)'; '*.*','All Files (*.*)'},'Select Deuterium Cell Image');

% Select corresponding spectral data 
uiwait(msgbox('Select the corresponding MATLAB spectral data to import for image analysis and processing', 'Select Spectral Data', 'modal'));
[Raman_Image_File, Raman_Image_Path] = uigetfile({'*.mat',...
    'MATLAB Files (*.mat)'; '*.*','All Files (*.*)'},'Select Spectral Data');


% -------------------------------------------------------------------------
%  2. Import, filter, and conduct blob analysis for cell image
% -------------------------------------------------------------------------
% Import cell image and filter
Cell_Image = imread(strcat(Cell_Image_Path, Cell_Image_File));
Cell_Image_Gray = rgb2gray(Cell_Image);
Cell_Image_BW = im2bw(Cell_Image_Gray,0.05);
Cell_Image_BW = bwareaopen(Cell_Image_BW,10);
Cell_Image_BW = imfill(Cell_Image_BW, 'holes');

% Create black border to remove image artefacts 
Cell_Image_BW(1:size(Cell_Image_BW,1), 1) = 0;
Cell_Image_BW(1:size(Cell_Image_BW,1), size(Cell_Image_BW,2)) = 0;
Cell_Image_BW(1, 1:size(Cell_Image_BW,1)) = 0;
Cell_Image_BW(size(Cell_Image_BW,1), 1:size(Cell_Image_BW,2)) = 0;

% Remove additional image artefacts
Cell_Row_Sums = sum(Cell_Image_BW,2);
Filter_Rows = Cell_Row_Sums > 0.95*size(Cell_Image_BW,2);
Cell_Column_Sums = sum(Cell_Image_BW,1);
Filter_Columns = Cell_Column_Sums > 0.95*size(Cell_Image_BW,1);
Cell_Image_BW(Filter_Rows,:) = 0;
Cell_Image_BW(:,Filter_Columns) = 0;

% Conduct blob analysis
Cell_Image_BW = bwareafilt(Cell_Image_BW,1);
Cell_Image_Perimeter = bwperim(Cell_Image_BW);
Cell_Image_Stats = regionprops(Cell_Image_BW,'Centroid','Area', 'PixelList');


% -------------------------------------------------------------------------
%  3. Import, filter, and conduct blob analysis for deuterium image
% -------------------------------------------------------------------------
% Import deuterium image and filter
Deuterium_Image = imread(strcat(Deuterium_Image_Path, Deuterium_Image_File));
Deuterium_Image_Gray = rgb2gray(Deuterium_Image);
Deuterium_Image_BW = im2bw(Deuterium_Image_Gray,0.2);
Deuterium_Image_BW = bwareaopen(Deuterium_Image_BW,10);

% Create black border to remove image artefacts 
Deuterium_Image_BW(1:size(Deuterium_Image_BW,1), 1) = 0;
Deuterium_Image_BW(1:size(Deuterium_Image_BW,1), size(Deuterium_Image_BW,2)) = 0;
Deuterium_Image_BW(1, 1:size(Deuterium_Image_BW,1)) = 0;
Deuterium_Image_BW(size(Deuterium_Image_BW,1), 1:size(Deuterium_Image_BW,2)) = 0;

% Remove additional image artefacts
Deuterium_Row_Sums = sum(Deuterium_Image_BW,2);
Filter_Rows = Deuterium_Row_Sums > 0.95*size(Deuterium_Image_BW,2);
Deuterium_Column_Sums = sum(Deuterium_Image_BW,1);
Filter_Columns = Deuterium_Column_Sums > 0.95*size(Deuterium_Image_BW,1);
Deuterium_Image_BW(Filter_Rows,:) = 0;
Deuterium_Image_BW(:,Filter_Columns) = 0;

% Conduct blob analysis
Deuterium_Image_Stats = regionprops(Deuterium_Image_BW,'Centroid','Area', 'PixelList');


% -------------------------------------------------------------------------
%  4. Import Raman spectral data
% -------------------------------------------------------------------------
% Import Raman spectral data
Raman_Image_Dataset = load(strcat(Raman_Image_Path, Raman_Image_File));
fn = fieldnames(Raman_Image_Dataset);
Raman_Image = reshape(double(Raman_Image_Dataset.(fn{1}).data), Raman_Image_Dataset.(fn{1}).imagesize(1), Raman_Image_Dataset.(fn{1}).imagesize(2), 1383);


% -------------------------------------------------------------------------
%  5. Generate deuterium-cell overlay image and get deuterium coordinates
% -------------------------------------------------------------------------
% Display deuterium-cell overlay image
Overlay_Image = cat(3, 255* uint8(Deuterium_Image_BW), 255* uint8(Cell_Image_BW), 0* uint8(Cell_Image_BW));

% Find deuterium pixels inside and outside the cell region
Deuterium_Pixel_List = [];
for i = 1:size(Deuterium_Image_Stats,1)
    Deuterium_Pixel_List = [Deuterium_Pixel_List; Deuterium_Image_Stats(i).PixelList];
end

% Get interior deuterium coordinates
Interior_Deuterium_Indices = ismember(Deuterium_Pixel_List, Cell_Image_Stats.PixelList, 'rows');
Interior_Deuterium_Coordinates = Deuterium_Pixel_List(Interior_Deuterium_Indices, :);

% Get exterior deuterium coordinates
Exterior_Deuterium_Indices = ~ismember(Deuterium_Pixel_List, Cell_Image_Stats.PixelList, 'rows');
Exterior_Deuterium_Coordinates = Deuterium_Pixel_List(Exterior_Deuterium_Indices, :);


% -------------------------------------------------------------------------
%  6. Calculate deuterium internalisation percentage and distances
% -------------------------------------------------------------------------
% Calculate interior and exterior deuterium percentage
Interior_Deuterium_Percentage = size(Interior_Deuterium_Coordinates,1) / size(Deuterium_Pixel_List,1) * 100;
Exterior_Deuterium_Percentage = size(Exterior_Deuterium_Coordinates,1) / size(Deuterium_Pixel_List,1) * 100;

% Calculate deuterium cell percentage
Deuterium_Cell_Percentage = size(Interior_Deuterium_Coordinates,1) / size(Cell_Image_Stats.PixelList,1) * 100;
Non_Deuterium_Cell_Percentage = 100 - Deuterium_Cell_Percentage;

% Calculate distance of each interior deuterium pixel from nearest cell
% perimeter pixel
[Cell_Perimeter_X, Cell_Perimeter_Y] = find(Cell_Image_Perimeter);
Cell_Perimeter_Coordinates = [Cell_Perimeter_Y, Cell_Perimeter_X];
[Deuterium_Distance_Indices, Deuterium_Distances] = dsearchn(Cell_Perimeter_Coordinates, Interior_Deuterium_Coordinates);


% -------------------------------------------------------------------------
%  7. Get Raman spectra for each deuterium pixel
% -------------------------------------------------------------------------
% Get all deuterium Raman spectra from Witec ProjectFOUR image matrix
Deuterium_Spectra_Coordinates = unique(round(Deuterium_Pixel_List * (size(Raman_Image,1) / size(Cell_Image,1))), 'rows');

Deuterium_Spectra = zeros(length(Deuterium_Spectra_Coordinates), 1383);

for i  = 1:size(Deuterium_Spectra_Coordinates,1)
    if Deuterium_Spectra_Coordinates(i,1) == 0
        Deuterium_Spectra_Coordinates(i,1) = 1;
    end
    if Deuterium_Spectra_Coordinates(i,2) == 0
        Deuterium_Spectra_Coordinates(i,2) = 1;
    end
    
    if Deuterium_Spectra_Coordinates(i,1) >= size(Raman_Image,1)
        Deuterium_Spectra_Coordinates(i,1) = size(Raman_Image,1);
    end
    if Deuterium_Spectra_Coordinates(i,2) >= size(Raman_Image,2)
        Deuterium_Spectra_Coordinates(i,2) = size(Raman_Image,2);
    end
        Deuterium_Spectra(i,:) = squeeze(Raman_Image(Deuterium_Spectra_Coordinates(i,2), Deuterium_Spectra_Coordinates(i,1),:));
end

% Get interior deuterium Raman spectra from Witec ProjectFOUR image matrix
Deuterium_Interior_Spectra_Coordinates = unique(round(Interior_Deuterium_Coordinates *(size(Raman_Image,1) / size(Cell_Image,1))), 'rows');

Deuterium_Interior_Spectra = zeros(length(Deuterium_Interior_Spectra_Coordinates), 1383);

for i  = 1:size(Deuterium_Interior_Spectra_Coordinates,1)
    if Deuterium_Interior_Spectra_Coordinates(i,1) == 0
        Deuterium_Interior_Spectra_Coordinates(i,1) = 1;
    end
    if Deuterium_Interior_Spectra_Coordinates(i,2) == 0
        Deuterium_Interior_Spectra_Coordinates(i,2) = 1;
    end
    
    if Deuterium_Interior_Spectra_Coordinates(i,1) >= size(Raman_Image,1)
        Deuterium_Interior_Spectra_Coordinates(i,1) = size(Raman_Image,1);
    end
    if Deuterium_Interior_Spectra_Coordinates(i,2) >= size(Raman_Image,2)
        Deuterium_Interior_Spectra_Coordinates(i,2) = size(Raman_Image,2);
    end
        Deuterium_Interior_Spectra(i,:) = squeeze(Raman_Image(Deuterium_Interior_Spectra_Coordinates(i,2), Deuterium_Interior_Spectra_Coordinates(i,1),:));
end


% Get exterior deuterium Raman spectra from Witec ProjectFOUR image matrix
Deuterium_Exterior_Spectra_Coordinates = unique(round(Exterior_Deuterium_Coordinates *(size(Raman_Image,1) / size(Cell_Image,1))), 'rows');

Deuterium_Exterior_Spectra = zeros(length(Deuterium_Exterior_Spectra_Coordinates), 1383);

for i  = 1:size(Deuterium_Exterior_Spectra_Coordinates,1)
    if Deuterium_Exterior_Spectra_Coordinates(i,1) == 0
        Deuterium_Exterior_Spectra_Coordinates(i,1) = 1;
    end
    if Deuterium_Exterior_Spectra_Coordinates(i,2) == 0
        Deuterium_Exterior_Spectra_Coordinates(i,2) = 1;
    end
    
    if Deuterium_Exterior_Spectra_Coordinates(i,1) >= size(Raman_Image,1)
        Deuterium_Exterior_Spectra_Coordinates(i,1) = size(Raman_Image,1);
    end
    if Deuterium_Exterior_Spectra_Coordinates(i,2) >= size(Raman_Image,2)
        Deuterium_Exterior_Spectra_Coordinates(i,2) = size(Raman_Image,2);
    end
        Deuterium_Exterior_Spectra(i,:) = squeeze(Raman_Image(Deuterium_Exterior_Spectra_Coordinates(i,2), Deuterium_Exterior_Spectra_Coordinates(i,1),:));
end

% Import X-axis for spectral plotting
load('C:\Users\ch215\Desktop\Data\axis_532.mat');


% -------------------------------------------------------------------------
%  8. Generate images for result visualisation
% -------------------------------------------------------------------------
% Plot segmented cell overlay image
figure(1)
imshow(Overlay_Image)
hold on
plot(Cell_Perimeter_Y, Cell_Perimeter_X, 'w.')
hold on
plot(Exterior_Deuterium_Coordinates(:,1), Exterior_Deuterium_Coordinates(:,2), 'b.')
hold on
plot(Interior_Deuterium_Coordinates(:,1), Interior_Deuterium_Coordinates(:,2), 'r.')
hold on
plot(Cell_Perimeter_Coordinates(Deuterium_Distance_Indices,1), Cell_Perimeter_Coordinates(Deuterium_Distance_Indices,2), 'm.')
legend('Cell Perimeter','Exterior Deuterium Signal','Interior Deuterium Signal', 'Nearest Cell Perimeter to Interior Deuterium Signal')

% Show image processing stages and EV analysis results
figure(2)
subplot(3,3,1), imshow(Cell_Image), title('Cell Image');
subplot(3,3,2), imshow(Cell_Image_Gray), title('Cell Image (Grayscale)');
subplot(3,3,3), imshow(Cell_Image_BW), title('Cell Image (Binary)');
subplot(3,3,4), imshow(Deuterium_Image), title('Deuterium Image');
subplot(3,3,5), imshow(Deuterium_Image_Gray), title('Deuterium Image (Grayscale)');
subplot(3,3,6), imshow(Deuterium_Image_BW), title('Deuterium Image (Binary)');
subplot(3,3,7), pie([Interior_Deuterium_Percentage, Exterior_Deuterium_Percentage]), title('EV Distribution'), legend('Interalised', 'External');
subplot(3,3,8), histogram(Deuterium_Distances), title('EV Internalisation Distances'), xlabel('Distance (pixels)'), ylabel('Counts');
subplot(3,3,9), pie([Deuterium_Cell_Percentage, Non_Deuterium_Cell_Percentage]), title('Normalised EV Internalisation'), legend('Deuterium Cell Content', 'Non-Deuterium Cell Content');

% Plot deuterium Raman spectra
figure(3)
subplot(2,1,1), plot(axis_532, mean(Deuterium_Spectra), axis_532, mean(Deuterium_Interior_Spectra), axis_532, mean(Deuterium_Exterior_Spectra)), xlim([500 3750]), title('Entire Spectrum')
subplot(2,1,2), plot(axis_532(1:850), mean(Deuterium_Spectra(:,1:850)), axis_532(1:850), mean(Deuterium_Interior_Spectra(:,1:850)), axis_532(1:850), mean(Deuterium_Exterior_Spectra(:,1:850))), xlim([500 2500]), title('Fingerprint and Silent Region'), legend('Mean of All Deuterium Spectra', 'Mean of Interior Deuterium Spectra', 'Mean of Exterior Deuterium Spectra')


% -------------------------------------------------------------------------
%  9. Save EV analysis information into struct array for further use
% -------------------------------------------------------------------------
Cell_EV_Analysis_Information.Interior_Deuterium_Percentage = Interior_Deuterium_Percentage;
Cell_EV_Analysis_Information.Exterior_Deuterium_Percentage = Exterior_Deuterium_Percentage;
Cell_EV_Analysis_Information.Deuterium_Distances = Deuterium_Distances;
Cell_EV_Analysis_Information.Deuterium_Cell_Percentage = Deuterium_Cell_Percentage;
Cell_EV_Analysis_Information.Non_Deuterium_Cell_Percentage = Non_Deuterium_Cell_Percentage;
Cell_EV_Analysis_Information.Deuterium_Interior_Spectra = Deuterium_Interior_Spectra;
Cell_EV_Analysis_Information.Deuterium_Exterior_Spectra = Deuterium_Exterior_Spectra;



