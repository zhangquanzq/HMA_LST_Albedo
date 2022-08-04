%% HMA Glcaier TimeSeries

%% Flags and preseted parameters.
% Flag assign the Data Type. 1 means MOD10A1, 2 means MYD10A1, 3 means MOD11A1, 4 means MYD11A1.
flg1 = 4;

% Flag assgin Day or night. 1 means Daytime, 2 means Nighttime.
flg2 = 1;

dataType = {'MOD10A1', 'MYD10A1', 'MOD11A1', 'MYD11A1'};
extentName = {'HMA', 'HMA', 'CN', 'CN'};
daynight = {'Day', 'Night'};

yearList = 2002:2019;
yearListN = length(yearList);

%% Paths.
% Root Paths.
dataRootPath = 'F:\HMA_LST_Albedo\Data';
dataTypeFolder = dataType{flg1};
dataStep2Folder = [dataType{flg1}, '_2_Mosaic', extentName{flg1}, '_TIF'];
dataStep3Folder = [dataType{flg1}, '_3_GlacierHMA_TIF'];
dataStep2Path = fullfile(dataRootPath, dataTypeFolder, dataStep2Folder);
dataStep3Path = fullfile(dataRootPath, dataTypeFolder, dataStep3Folder);

dataYearFolderList = dir(fullfile(dataStep2Path, '*_TIF'));
dataYearFolderList = {dataYearFolderList.name}';
dataYearFolderListN = length(dataYearFolderList);

% Get the dataYearFolder to be handled based on yearList.
dataYearList = zeros(dataYearFolderListN, 1);
for i = 1 : dataYearFolderListN
    dataYearFolder = dataYearFolderList{i};
    dataYearList(i) = str2double(dataYearFolder(9:12));
end
index = zeros(yearListN, 1);
for i = 1 : yearListN
    index(i) = find(dataYearList == yearList(i));
end
dataYearFolderList = dataYearFolderList(index);
dataYearFolderListN = length(dataYearFolderList);

%% Read references and attributes.
% Read the spatial reference of HMA_Glacier_Extent images.
glacierExtentPath = fullfile(dataRootPath, 'Basic\Raster\HMA_Glacier_Extent.tif');
[glacierLayer, glacierReference] = geotiffread(glacierExtentPath);
glacierLayer = double(glacierLayer);
glacierLayer(glacierLayer == 0) = nan;
glacierRowN = glacierReference.RasterSize(1);
glacierColN = glacierReference.RasterSize(2);
glacierCellsizeX = glacierReference.CellExtentInWorldX;
glacierCellsizeY = glacierReference.CellExtentInWorldY;
glacierXMin = glacierReference.XWorldLimits(1);
glacierXMax = glacierReference.XWorldLimits(2);
glacierYMin = glacierReference.YWorldLimits(1);
glacierYMax = glacierReference.YWorldLimits(2);
glacierPixelLeftSideXVector = glacierXMin : glacierCellsizeX : glacierXMax - glacierCellsizeX;
glacierPixelTopSideYVector = glacierYMax : - glacierCellsizeY : glacierYMin + glacierCellsizeY;

% Read the spatial reference of MODIS LST images.

modisPath = fullfile(dataStep2Path, [dataType{flg1}, ['_2008XXX_TIF\', dataType{flg1}, '.A2008001.LST_Day.tif']]);
[~, modisReference] = geotiffread(modisPath);
modisRowN = modisReference.RasterSize(1);
modisColN = modisReference.RasterSize(2);
modisCellsizeX = modisReference.CellExtentInWorldX;
modisCellsizeY = modisReference.CellExtentInWorldY;
modisXMin = modisReference.XWorldLimits(1);
modisXMax = modisReference.XWorldLimits(2);
modisYMin = modisReference.YWorldLimits(1);
modisYMax = modisReference.YWorldLimits(2);
modisPixelLeftSideXVector = modisXMin : modisCellsizeX : modisXMax - modisCellsizeX;
modisPixelTopSideYVector = modisYMax : - modisCellsizeY : modisYMin + modisCellsizeY;

% Get the starting/ending row numbers and column numbers from the modis images that are aligned with
%   the boundary of HMA_Glacier_Extent images.
if glacierXMin > modisXMin
    modisStartCol = ceil((glacierXMin - modisXMin) / modisCellsizeX);
else
    modisStartCol = 1;
end

if modisYMax > glacierYMax
    modisStartRow = ceil((modisYMax - glacierYMax) / modisCellsizeY);
else
    modisStartRow = 1;
end

if glacierXMax < modisXMax
    modisEndCol = modisColN - ceil((modisXMax - glacierXMax) / modisCellsizeX) + 1;
else
    modisEndCol = modisColN;
end

if modisYMin < glacierYMin
    modisEndRow = modisRowN - ceil((glacierYMin - modisYMin) / modisCellsizeY) + 1;
else
    modisEndRow = modisRowN;
end

%% Operate the files in each year folder.
% dateN = yearListN * 365 + sum(leapyear(yearList));
% [dateList, hmaGlacierMeanLstList] = deal(zeros(dateN, 1));
[dateList, hmaGlacierLstMeanList, hmaGlacierLstStdList] = deal([]);
for i = 1 : dataYearFolderListN
    % Read the LST and QC file list.
    dataYearPath = fullfile(dataStep2Path, dataYearFolderList{i});
    
    lstYearDateName = [dataType{flg1}, '*LST_', daynight{flg2}, '.tif'];
    lstYearDateList = dir(fullfile(dataYearPath, lstYearDateName));
    lstYearDateList = {lstYearDateList.name}';
    lstYearDateListN = length(lstYearDateList);
    
    qcYearDateName = [dataType{flg1}, '*QC_', daynight{flg2}, '.tif'];
    qcYearDateList = dir(fullfile(dataYearPath, qcYearDateName));
    qcYearDateList = {qcYearDateList.name}';
    qcYearDateListN = length(qcYearDateList);
    
    if lstYearDateListN ~= qcYearDateListN
        disp(['Check the data completeness in Folder: ', dataYearFolderList]);
        break
    end
    
    for j = 1 : lstYearDateListN
        lstYearDateName = lstYearDateList{j};
        [lstLayer, lstRef] = geotiffread(fullfile(dataYearPath, lstYearDateList{j}));
        qcLayer = geotiffread(fullfile(dataYearPath, qcYearDateList{j}));
        lstRowN = lstRef.RasterSize(1);
        lstColN = lstRef.RasterSize(2);
        if lstRowN ~= modisRowN || lstColN ~= modisColN
            continue
        end
        lstLayer = double(lstLayer(modisStartRow:modisEndRow, modisStartCol:modisEndCol)) * 0.02;
        qcLayer = double(qcLayer(modisStartRow:modisEndRow, modisStartCol:modisEndCol));
        lstLayer(qcLayer ~= 0) = nan;
        lstGlacierLayer = lstLayer .* glacierLayer;
        dateList(end+1) = str2double(lstYearDateName(10:16));
        hmaGlacierLstMeanList(end+1) = mean(lstGlacierLayer(:), 'omitnan');
        hmaGlacierLstStdList(end+1) = std(lstGlacierLayer(:), 'omitnan');
        lstGlacierFolder = fullfile(dataStep3Path, dataYearFolderList{i});
        if ~exist(lstGlacierFolder, 'dir')
            mkdir(lstGlacierFolder);
        end
        glacierLstName = replace(lstYearDateName, 'LST', 'HMA_Glacier_LST');
        lstGlacierPath = fullfile(lstGlacierFolder, glacierLstName);
        if  ~exist(lstGlacierPath, 'file')
            geotiffwrite(lstGlacierPath, lstGlacierLayer, modisReference);
        end
        disp(lstYearDateName);
    end
end
lstArrayMatName = [dataType{flg1}, '_LST_sta.csv'];
lstArrayMatPath = fullfile(dataStep3FolderPath, lstArrayMatName);
lstStatisticsArray = [dateList', hmaGlacierLstMeanList', hmaGlacierLstStdList'];
writematrix(lstStatisticsArray, lstArrayMatPath);

dateList2 = {};
for i = 1: length(dateList)
    dateList2{end+1} = yday2ymd(num2str(dateList(i)));
end



















