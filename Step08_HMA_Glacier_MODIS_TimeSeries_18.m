%% HMA冰川的MODIS地表温度时间序列.

% !!! 坡度, 高程分级标准目前是等间距, 导致每个间隔内的像元样本数量差异大，看是否能改成保证每个间隔内数量一致
%   的分级标准. !!!

%% 标记和预设参数.
% 指定数据类型的标记. 1表示MOD10A1, 2表示MYD10A1, 3表示MOD11A1, 4表示MYD11A1.
flg1 = 4;
% 指定昼夜的标记. 1表示白天, 2表示晚上.
flg2 = 1;
% 指定面积比例的标记. 1表示100, 2表示95, 3表示90, 4表示85, 5表示80.
flg3 = 5;
% 指定时间序列图横轴坐标类型的标记. 1表示月, 2表示日.
flg4 = 1;

% Region list = {'Altay and Sayan', 'C Himalaya', 'E Himalaya', 'E Kun Lun (Altyn Tagh)',...
%     'E Tien Shan (Dzhungaria)', 'Hengduan Shan', 'Hindu Kush', 'Hissar Alay', 'Inner Tibet',...
%     'Karakoram', 'Pamir (Safed Khirs/W Tarim)', 'Qilian Shan', 'S and E Tibet', 'W Himalaya',...
%     'W Kun Lun', 'W Tien Shan'};
% 季节划分 (冬: 12,1,2, 春: 3,4,5, 夏: 6,7,8, 秋: 9,10,11)

% 数据类型, 名称, 范围, 分辨率, 昼夜, 统计时间间隔.
dataType = {'MOD10A1', 'MYD10A1', 'MOD11A1', 'MYD11A1'};
dataType = dataType{flg1};

dataName = {'Albedo', 'LST'};
dataName = dataName{round(flg1/2)};

extent = {'HMA', 'CN'};
extent = extent{round(flg1/2)};

cellsize = {'500m', '1km'};
cellsize = cellsize{round(flg1/2)};

daynight = {'Day', 'Night'};
daynight = daynight{flg2};

dateType = {'Month', 'Day'};
dateType = dateType{flg4};

% 年份列表, 季节, 像元面积百分比.
yearList = {2002 : 2019, 2000 : 2019};
yearList = yearList{mod(flg1, 2) + 1};
yearListN = length(yearList);

seasons = struct('Winter1', [1, 2], 'Spring', [3, 4, 5], 'Summer', [6, 7, 8], 'Autumn', ...
    [9, 10, 11], 'Winter2', 12);

pctList = {'100', '95', '90', '85', '80'};
pct = pctList{flg3}; minPct = pctList{end};


% 高程和坡度区间.
elevMeanEdges = 2000: 200: 8600;
slpMeanEdges = 0: 10: 90;

% LST的有效qc取值.
qcValid = [0, 17];

%% 路径.
% 根目录.
rootDir = 'G:\HMA_LST_Albedo\Data';
stepsDir = fullfile(rootDir, 'GlacierAreaInPixel');
modisDir = fullfile(rootDir, dataType);

% 各输入数据存放的文件夹路径.
% hmaPixelStaDir = fullfile(stepsDir, 'Step7_HMA_Pixel_Sta');
hmaPixelStaPctDir = fullfile(stepsDir, 'Step8_HMA_Pixel_Sta_Percent');
hmaPixelPctRasterDir = fullfile(stepsDir, 'Step9_HMA_Pixel_Percent_Raster');

% 输出Mat文件的文件夹路径.
hmaMatDir = fullfile(stepsDir, 'Step10_HMA_Matlab');
if ~exist(hmaMatDir, 'dir')
    mkdir(hmaMatDir)
end

tableDir = fullfile(stepsDir, 'Step11_HMA_Tables');
figDir = fullfile(stepsDir, 'Step12_HMA_Figures');

% 最低冰川面积比例(80%)以上的像元FID编号的栅格文件.
hmaMinPctRasterName = ['HMA_ModisPixel_', cellsize, '_rgi60_FID_', minPct, 'percent.tif'];
hmaMinPctRasterPath = fullfile(hmaPixelPctRasterDir, hmaMinPctRasterName);

% 各冰川面积比例像元FID编号的栅格文件.
hmaPctRasterName = ['HMA_ModisPixel_', cellsize, '_rgi60_FID_', pct, 'percent.tif'];
hmaPctRasterPath = fullfile(hmaPixelPctRasterDir, hmaPctRasterName);

% 每个像元内高程, 坡度统计的矢量文件.
hmaMinPctStaElevName = ['HMA_ModisPixel_', cellsize, '_sta_elev_', minPct, 'percent.shp'];
hmaMinPctStaElevPath = fullfile(hmaPixelStaPctDir, hmaMinPctStaElevName);

hmaMinPctStaSlpName = ['HMA_ModisPixel_', cellsize, '_sta_slp_', minPct, 'percent.shp'];
hmaMinPctStaSlpPath = fullfile(hmaPixelStaPctDir, hmaMinPctStaSlpName);

% hmaPixelStaElevName = ['HMA_ModisPixel_', cellsize, '_sta_elev_', pct, 'percent.shp'];
% hmaPixelStaElevPath = fullfile(hmaPixelStaPctDir, hmaPixelStaElevName);

% hmaPixelStaSlpName = ['HMA_ModisPixel_', cellsize, '_sta_slp_', pct, 'percent.shp'];
% hmaPixelStaSlpPath = fullfile(hmaPixelStaPctDir, hmaPixelStaSlpName);

% 每年的MODIS冰川像元值Mat文件名(用在sprintf函数中).
modisMatName = ['HMA_', dataType, '_%s_', minPct, 'percent.mat'];

% 获取最低冰川面积比例(80%)以上像元的位置索引和FID值.
[hmaMinPctLayer, hmaMinPctRef] = readgeoraster(hmaMinPctRasterPath);
minPctNodataValue = georasterinfo(hmaMinPctRasterPath).MissingDataIndicator;
geoTag = geotiffinfo(hmaMinPctRasterPath).GeoTIFFTags.GeoKeyDirectoryTag;

[hmaPctLayer, hmaPctRef] = readgeoraster(hmaPctRasterPath);
pctNodataValue = georasterinfo(hmaPctRasterPath).MissingDataIndicator;

minPctIndexLayer = (hmaMinPctLayer ~= minPctNodataValue);
minPctPixelN = sum(minPctIndexLayer(:));
minPctFidList = hmaMinPctLayer(minPctIndexLayer);
[hmaRowN, hmaColN] = size(hmaMinPctLayer);

pctIndexLayer = (hmaPctLayer ~= pctNodataValue);
pctPixelN = sum(pctIndexLayer(:));
pctFidList = hmaPctLayer(pctIndexLayer);

% 获取冰川MODIS像元的高程, 坡度统计表格, 按MODIS Pixel FID对属性表排序.
[~, staElevTable] = shaperead(hmaMinPctStaElevPath);
staElevTable = staElevTable(minPctFidList + 1);  % FID从0开始，加1从1开始.
elevMeanRecords = [staElevTable.MEAN]';
elevPctRecords = [staElevTable.Area_Pct]';

[~, staSlpTable] = shaperead(hmaMinPctStaSlpPath);
staSlpTable = staSlpTable(minPctFidList + 1);  % FID从0开始，加1从1开始.
slpMeanRecords = [staSlpTable.MEAN]';
slpPctRecords = [staSlpTable.Area_Pct]';

% 获取冰川分区名称列表.
regionRecords = {staElevTable.FULL_NAME}';
regionList = unique(regionRecords);
regionListN = length(regionList);

% HMA冰川区LST或Albedo栅格文件目录.
modisMosaicFolder = sprintf('%s_2_Mosaic%s_TIF', dataType, extent);

modisPctFolder = sprintf('%s_3_HMAGlacier%s_TIF', dataType, pct);

modisMeanDir = fullfile(modisDir, sprintf('%s_4_HMAGlacier%s_Mean_TIF', dataType, pct));
if ~exist(modisMeanDir, 'dir')
    mkdir(modisMeanDir);
end

%% 将每一年HMA冰川区的MODIS数据存储为Mat文件, 提高程序运行效率.
for i = 1 : yearListN
    break
    yearStr = num2str(yearList(i));
    
    % 判断Mat文件是否存在.
    fmtStr = {yearStr, [daynight, '_', yearStr]};
    modisMatPath = fullfile(hmaMatDir, sprintf(modisMatName, fmtStr{round(flg1/2)}));
    if exist(modisMatPath, 'file')
        continue
    end

    % 读取MODIS LST/Albedo数据的文件名列表.
    modisYearFolder = [dataType, '_', yearStr, 'XXX_TIF'];
    modisTifName = {[dataType, '*Albedo.tif'], [dataType, '*LST_', daynight, '.tif']};
    modisYearDir = fullfile(modisDir, modisMosaicFolder, modisYearFolder);
    modisTifList = dir(fullfile(modisYearDir, modisTifName{round(flg1/2)}));
    modisTifList = {modisTifList.name}';
    modisTifListN = length(modisTifList);
    
    % 若处理LST数据, 则再读取QC文件名列表, 并检查LST和QC文件是否对应.
    if strcmp(dataName, 'LST')
        qcTifName = [dataType, '*QC_', daynight, '.tif'];
        qcTifList = dir(fullfile(modisYearDir, qcTifName));
        qcTifList = {qcTifList.name}';
        if modisTifListN ~= length(qcTifList)
            error(['LST和QC文件数量不一致，请检查文件夹：', modisYearFolder]);
        end
    end
    
    % 读取MODIS TIF格式文件, 并从数据文件名中获取数据日期列表, 存到对应变量中, modisMatrix和qcMatrix为2D
    %   矩阵, 行为像元, 列为日期.
    modisMatrix = zeros(minPctPixelN, modisTifListN, 'single') * nan;
    if strcmp(dataName, 'LST')
        qcMatrix = ones(minPctPixelN, modisTifListN, 'uint16') * 65535;
    end
    modisDateList = cell(modisTifListN, 1);
    for j = 1 : modisTifListN
        modisTifName = modisTifList{j};
        modisDateList{j} = yday2ymd(modisTifName(10:16));
        modisTifPath = fullfile(modisYearDir, modisTifName);
        modisLayer = single(readgeoraster(modisTifPath));
        modisNodata = single(georasterinfo(modisTifPath).MissingDataIndicator);

        if strcmp(dataName, 'LST') % 处理LST数据.
            qcLayer = readgeoraster(fullfile(modisYearDir, qcTifList{j}));
            if isequal(size(modisLayer), size(qcLayer), [hmaRowN, hmaColN])
                modisLayer(ismember(modisLayer, [modisNodata, 0, 65535])) = nan;
                modisMatrix(:, j) = modisLayer(minPctIndexLayer) * 0.02;
                qcMatrix(:, j) = qcLayer(minPctIndexLayer);
                disp(['添加图幅到Mat文件: ', modisTifName]);
            else
                disp(['不完整图幅: ', modisTifName]);
            end
        elseif strcmp(dataName, 'Albedo') % 处理Albedo数据.
            if isequal(size(modisLayer), [hmaRowN, hmaColN])
                modisLayer(modisLayer == modisNodata) = nan;
                modisMatrix(:, j) = modisLayer(minPctIndexLayer);
                disp(['添加图幅到Mat文件: ', modisTifName]);
            else
                disp(['不完整图幅: ', modisTifName]);
            end
        end
    end
    modisDateList = datetime(modisDateList, 'InputFormat', 'yyyy/MM/dd');

    % 保存变量到Mat文件, 存储了HMA冰川区每年的LST/Albedo数据, qc数据, 以及时间列表.
    save(modisMatPath, 'modisMatrix', 'modisDateList');
    if strcmp(dataName, 'LST')
        save(modisMatPath, 'qcMatrix', '-append');
    end
end

%% 计算并输出MODIS冰川区LST/Albedo数据的年, 月平均值影像.
for i = 1 : yearListN
    yearStr = num2str(yearList(i));
    
    % 读取Mat文件, 并对LST或Albedo做质量控制.
    fmtStr = {yearStr, [daynight, '_', yearStr]};
    modisMatPath = fullfile(hmaMatDir, sprintf(modisMatName, fmtStr{round(flg1/2)}));
    load(modisMatPath, 'modisMatrix', 'modisDateList');  
    if strcmp(dataName, 'LST')
        load(modisMatPath, 'qcMatrix');
        modisMatrix(~ismember(qcMatrix, qcValid)) = nan;
    elseif strcmp(dataName, 'Albedo')
        modisMatrix(modisMatrix == 0 | modisMatrix > 100) = nan;
    end

    % 创建输出TIF文件的年文件夹.
    modisYearFolder = [dataType, '_', yearStr, 'XXX_TIF'];
    modisYearDir = fullfile(modisMeanDir, modisYearFolder);
    if ~exist(modisYearDir, 'dir')
        mkdir(modisYearDir)
    end

    % MODIS冰川区LST/Albedo的年平均值.
    modisYearMeanName = [dataType, '.A', yearStr, '.', dataName, '_', daynight, '.tif'];
    modisYearMeanPath = fullfile(modisYearDir, modisYearMeanName);
    if ~exist(modisYearMeanPath, 'file')
        modisYearMeanLayer = zeros(hmaRowN, hmaColN, 'single');
        modisYearMeanLayer(minPctIndexLayer) = mean(modisMatrix, 2, 'omitnan');
        modisYearMeanLayer(~pctIndexLayer) = nan;
        geotiffwrite(modisYearMeanPath, modisYearMeanLayer, hmaMinPctRef, ...
            GeoKeyDirectoryTag=geoTag, tifftags=struct('Compression','LZW'));
        disp(['输出: ', modisYearMeanName]);
    end
    
    % MODIS冰川区LST/Albedo的月平均值.
    monthList = modisDateList.Month;
    monthType = unique(monthList);
    modisMonthMeanLayer = zeros(hmaRowN, hmaColN, 'single');
    for j = 1 : length(monthType)
        monthNum = monthType(j); monthStr = num2str(monthNum, '%02d');
        modisMonthMeanName = replace(modisYearMeanName, yearStr, [yearStr, monthStr]);
        modisMonthMeanPath = fullfile(modisYearDir, modisMonthMeanName);
        if exist(modisMonthMeanPath, 'file')
            continue
        end
        modisMonthMeanVector = mean(modisMatrix(:, monthList == monthNum), 2, 'omitnan');
        modisMonthMeanLayer(minPctIndexLayer) = modisMonthMeanVector;
        modisMonthMeanLayer(~pctIndexLayer) = nan;
        geotiffwrite(modisMonthMeanPath, modisMonthMeanLayer, hmaMinPctRef, ...
            GeoKeyDirectoryTag=geoTag, tifftags=struct('Compression','LZW'));
        disp(['输出: ', modisMonthMeanName]);
    end
end

%% 导出表格, 作图.
% 创建存储各时间序列表格, 图的文件夹.

dataPctFolder = [dataType, '_GlacierAreaPercent', pct];
xlsDir = fullfile(tableDir, dataPctFolder);
if ~exist(xlsDir, 'dir')
    mkdir(xlsDir);
end
figPctDir = fullfile(figDir, dataPctFolder);
if ~exist(figPctDir, 'dir')
    mkdir(figPctDir);
end

% 按年周期创建不同分区冰川表面温度的时间序列图, 按高程, 坡度统计.
hmaModisYearMeanVector = zeros(yearListN, 1);
hmaParaSeasonMeanMatrix = zeros(yearListN, 4);
for i = 1 : yearListN
    yearStr = num2str(yearList(i));
    
    % modisMat文件存储了HMA冰川区每年的LST/Albedo数据.
    fmtStr = {yearStr, [daynight, '_', yearStr]};
    modisMatPath = fullfile(hmaMatDir, sprintf(modisMatName, fmtStr{round(flg1/2)}));
    load(modisMatPath, 'modisMatrix', 'modisDateList');  % modisMatrix为2D矩阵, 行代表像元, 列为日期.
    
    hmaModisYearMeanVector(i) = mean(modisMatrix, 'all', 'omitnan');

    % 分区统计.
    for j =  1 : regionListN
        regionName = regionList{j};
        regionName2 = replace(regionName, '/', '_');

        % 将分级统计的高程结果存储为Excel文件。
        elevLstXlsName = sprintf('%stime_GST_timeseries_%s_%s_%s_Elev.xlsx', daynight, dataType, ...
            dateType, regionName2);
        elevAlbedoXlsName = sprintf('Albedo_timeseries_%s_%s_%s_Elev.xlsx', dataType, dateType, ...
            regionName2);
        elevXlsName = {elevAlbedoXlsName, elevLstXlsName};
        elevXlsName = elevXlsName{round(flg1/2)};
        elevXlsPath = fullfile(xlsDir, elevXlsName);
        elevXlsExist = exist(elevXlsPath, 'file');
        if ~elevXlsExist || (elevXlsExist && ~ismember(yearStr, sheetnames(elevXlsPath)))
            % 按高程分级统计.
            regionIndex = strcmp(regionName, regionRecords);
            subParaMatrix = modisMatrix(regionIndex, :);
            subElevMeanRecords = elevMeanRecords(regionIndex);
            [elevBinsRange, elevBinsCount, elevDateTypes, elevDateCount, paraArrayElevMean] = ...
                intervalMean(subParaMatrix, subElevMeanRecords, elevMeanEdges, modisDateList, ...
                dateType);
            elevBinsRangeN = length(elevBinsCount); elevDateTypeN = length(elevDateTypes);
            elevBinsRangeCell = cell(elevBinsRangeN, 1);
            for k = 1 : elevBinsRangeN
                elevBinsMin = elevBinsRange(k, 1); elevBinsMax = elevBinsRange(k, 2);
                elevBinsRangeCell(k) = {sprintf('%d-%d m', elevBinsMin, elevBinsMax)};
            end
            elevStaCell = cell(elevBinsRangeN + 2, elevDateTypeN + 2);
            elevStaCell(1, 2) = {dateType};
            elevStaCell(1, 3:end) = num2cell(elevDateTypes');
            elevStaCell(2, 1) = {'Elevation range'};
            elevStaCell(2, 2) = {'Count'};
            elevStaCell(2, 3:end) = num2cell(elevDateCount');
            elevStaCell(3:end, 1) = elevBinsRangeCell;
            elevStaCell(3:end, 2) = num2cell(elevBinsCount);
            elevStaCell(3:end, 3:end) = num2cell(paraArrayElevMean);
            writecell(elevStaCell, elevXlsPath, 'Sheet', yearStr);
            fprintf('添加%s年统计值到%s\n', yearStr, elevXlsName);
        end
        
        % 将分级统计的坡度结果存储为Excel文件.
        slpLstXlsName = sprintf('%stime_GST_timeseries_%s_%s_%s_Slp.xlsx', daynight, dataType, ...
            dateType, regionName2);
        slpAlbedoXlsName = sprintf('Albedo_timeseries_%s_%s_%s_Slp.xlsx', dataType, dateType, ...
            regionName2);
        slpXlsName = {slpAlbedoXlsName, slpLstXlsName};
        slpXlsName = slpXlsName{round(flg1/2)};
        slpXlsPath = fullfile(xlsDir, slpXlsName);
        slpXlsExist = exist(slpXlsPath, 'file');
        if ~slpXlsExist || (slpXlsExist && ~ismember(yearStr, sheetnames(slpXlsPath)))
            % 按坡度分级统计.
            regionIndex = strcmp(regionName, regionRecords);
            subParaMatrix = modisMatrix(regionIndex, :);
            subSlpMeanRecords = slpMeanRecords(regionIndex);
            [slpBinsRange, slpBinsCount, slpDateTypes, slpDateCount, paraArraySlpMean] = ...
                intervalMean(subParaMatrix, subSlpMeanRecords, slpMeanEdges, modisDateList, ...
                dateType);
            slpBinsRangeN = length(slpBinsCount); slpDateTypeN = length(slpDateTypes);
            slpBinsRangeCell = cell(slpBinsRangeN, 1);
            for k = 1 : slpBinsRangeN
                slpBinsMin = slpBinsRange(k, 1); slpBinsMax = slpBinsRange(k, 2);
                slpBinsRangeCell(k) = {sprintf('%d-%d °', slpBinsMin, slpBinsMax)};
            end
            slpStaCell = cell(slpBinsRangeN + 2, slpDateTypeN + 2);
            slpStaCell(1, 2) = {dateType};
            slpStaCell(1, 3:end) = num2cell(slpDateTypes');
            slpStaCell(2, 1) = {'Slope range'};
            slpStaCell(2, 2) = {'Count'};
            slpStaCell(2, 3:end) = num2cell(slpDateCount');
            slpStaCell(3:end, 1) = slpBinsRangeCell;
            slpStaCell(3:end, 2) = num2cell(slpBinsCount);
            slpStaCell(3:end, 3:end) = num2cell(paraArraySlpMean);
            writecell(slpStaCell, slpXlsPath, 'Sheet', yearStr);
            fprintf('添加%s年统计值到%s\n', yearStr, slpXlsName);
        end
        
        % 创建存储各时间序列图的文件夹.
        figRegionDir = fullfile(figPctDir, regionName2);
        if ~exist(figRegionDir, 'dir')
            mkdir(figRegionDir);
        end
        
        % 制图, 高程.
        f1LstName = sprintf('%stime_GST_timeseries_%s_%s_%s_%s_Elev.png', daynight, dataType, ...
            dateType, regionName2, yearStr);
        f1AlbedoName = sprintf('Albedo_timeseries_%s_%s_%s_%s_Elev.png', dataType, dateType, ...
            regionName2, yearStr);
        f1Name = {f1AlbedoName, f1LstName};
        f1Name = f1Name{round(flg1/2)};
        f1Path = fullfile(figRegionDir, f1Name);
        if ~exist(f1Path, 'file')
            f1 = gstFigure(elevDateTypes, paraArrayElevMean, elevBinsRange, elevBinsCount, ...
                dateType, daynight, dataType, regionName, yearStr, 'Elevation');
            exportgraphics(f1, f1Path);
%             print(f1, f1Path, '-dpng', '-r200');
        end
        
        % 制图，坡度。
        f2LstName = sprintf('%stime_GST_timeseries_%s_%s_%s_%s_Slp.png', daynight, dataType, ...
            dateType, regionName2, yearStr);
        f2AlbedoName = sprintf('Albedo_timeseries_%s_%s_%s_%s_Slp.png', dataType, dateType, ...
            regionName2, yearStr);
        f2Name = {f2AlbedoName, f2LstName};
        f2Name = f2Name{round(flg1/2)};
        f2Path = fullfile(figRegionDir, f2Name);
        if ~exist(f2Path, 'file')
            f2 = gstFigure(slpDateTypes, paraArraySlpMean, slpBinsRange, slpBinsCount, ...
                dateType, daynight, dataType, regionName, yearStr,  'Slope');
            exportgraphics(f2, f2Path);
%             print(f2, f2Path, '-dpng', '-r200');
        end
        close all;
    end
end

%% 季节和年纪变化趋势图.
% 季节变化趋势数据.
seasonParaValueMatrix = zeros(4, yearListN);
for i = 1 : yearListN
    yearStr1 = num2str(yearList(i));  % 当前年.
    fmtStr1 = {yearStr1, [daynight, '_', yearStr1]};
    paraMatPath1 = fullfile(hmaMatDir, sprintf(modisMatName, fmtStr1{round(flg1/2)}));
    paraMatStruct1 = load(paraMatPath1, 'paraMatrix', 'dateList');
    dateList1 = paraMatStruct1.dateList;
    paraMatrix1 = paraMatStruct1.paraMatrix;
    monthList1 = dateList1.Month;

    if i ~= 1
        yearStr2 = num2str(yearList(i-1));  % 上一年.
        fmtStr2 = {yearStr2, [daynight, '_', yearStr2]};
        paraMatPath2 = fullfile(hmaMatDir, sprintf(modisMatName, fmtStr2{round(flg1/2)}));
        paraMatStruct2 = load(paraMatPath2, 'paraMatrix', 'dateList');
        dateList2 = paraMatStruct2.dateList;
        paraMatrix2 = paraMatStruct2.paraMatrix;
        monthList2 = dateList2.Month;
    end

    if i == 1
        paraWinterMatrix = paraMatrix1(:, ismember(monthList1, seasons.Winter1));
    else
        paraWinter1Matrix = paraMatrix1(:, ismember(monthList1, seasons.Winter1));
        paraWinter2Matrix = paraMatrix2(:, ismember(monthList2, seasons.Winter2));
        paraWinterMatrix = [paraWinter1Matrix, paraWinter2Matrix];
    end
    paraSpringMatrix = paraMatrix1(:, ismember(monthList1, seasons.Spring));
    paraSummerMatrix = paraMatrix1(:, ismember(monthList1, seasons.Summer));
    paraAutumnMatrix = paraMatrix1(:, ismember(monthList1, seasons.Autumn));
    
    winterMean = mean(paraWinterMatrix, 'all', 'omitnan');
    sprintMean = mean(paraSpringMatrix, 'all', 'omitnan');
    summerMean = mean(paraSummerMatrix, 'all', 'omitnan');
    autumnMean = mean(paraAutumnMatrix, 'all', 'omitnan');
    seasonParaValueMatrix(:, i) = [winterMean, sprintMean, summerMean, autumnMean];
end

% 季节变化趋势图.
f3 = figure; f3.Position = [400 300 1000 400];
for i = 1 : 4
    seasonParaValueVector = seasonParaValueMatrix(i, :) - 273.15;
    p1 = polyfit(yearList, seasonParaValueVector, 1);
    plot(yearList, seasonParaValueVector, 'o-',  ...
        [2000, 2020], [p1(1) * 2000 + p1(2), p1(1) * 2020 + p1(2)], 'k');
    hold on;
end
legend('Winter', '', 'Spring', '', 'Summer', '', 'Autumn', '', 'Location', 'southeast');
titleStr = ['Seasonal Average of %s-time GST in HMA from %d to %d\n'...
    '(%s, Glacier pixel threshold: %s%%)'];
t = title(sprintf(titleStr, daynight, yearList(1), yearList(end), dataType, pct));
t.FontWeight = 'bold';
xlabel('Year'); ylabel('Glacier surface temperature (°C)');
f3NameStr = 'Seasonal Average of %s-time GST in HMA from %d to %d, %s, %s%%';
f3Name = sprintf(f3NameStr, daynight, yearList(1), yearList(end), dataType, pct);
print(f3, f3Name, '-dpng', '-r200');

% 年纪变化趋势图.
f3 = figure; f3.Position = [400 300 1000 400];
p = polyfit(yearList, hmaModisYearMeanVector - 273.15, 1);
plot(yearList, hmaModisYearMeanVector - 273.15, 'o-',  ...
    [2000, 2020], [p(1) * 2000 + p(2), p(1) * 2020 + p(2)], ...
    'MarkerFaceColor', 'b');
titleStr = ['Annual Average of %s-time GST in HMA from %d to %d\n'...
    '(%s, Glacier pixel threshold: %s%%)'];
t = title(sprintf(titleStr, daynight, yearList(1), yearList(end), dataType, pct));
t.FontWeight = 'bold';
xlabel('Year'); ylabel('Glacier surface temperature (°C)');
l = legend('Temperature', 'Trend');
l.Location = 'northwest';
f3NameStr = 'Annual Average of %s-time GST in HMA from %d to %d, %s, %s%%';
f3Name = sprintf(f3NameStr, daynight, yearList(1), yearList(end), dataType, pct);
print(f3, f3Name, '-dpng', '-r200');
close all;


%% 自定义函数。
% 按某一参数的统计值划分区间，计算每一区间的均值。
function [paramStaBinsRange, paramStaBinsCount, dateTypes,  dateIntervalCount,...
    lstArrayParamSta] = intervalMean(lstArray, paramStaRecords, paramStaIntervals, datetimeList,...
    dateIntervalFlag)
% 输入参数：
% lstArray 是冰川表面温度值矩阵，行表示不同的空间位置，列表示日期（年内）。
% paramStaRecords 是参数统计值列向量，每个值表示lstArray中每行对应的像元位置的参数统计。
% paramStaIntervals 是参数统计值的划分区间。例：坡度取值0-90度，可划分为5度间隔的多个区间（0:5:90）。
% datetimeList 是某年内的日期列表。
% dateIntervalFlag 是日期划分区间标记。可取'Month'和'Day'。

% 输出参数：
% paramStaBinsRange 是参数统计值区间范围矩阵，每行的2数字表示区间值域，行数为参数统计值区间的个数。
% paramStaBinsCount 是参数统计值区间内的统计值个数向量，行数为参数统计值区间的个数。
% dateIntervalTypes 是日期列表按dateIntervalFlag划分后的日期区间值。
% dateIntervalCount 是每个日期区间包含的日期个数。
% lstArrayParamSta 是按参数取值区间计算的地表温度值矩阵。
% 注意：paramStaRecords中相邻两个统计值的所属区间没有特定顺序。

paramStaBins = discretize(paramStaRecords, paramStaIntervals);
paramStaBinsType = unique(paramStaBins);
paramStaBinsTypeN = length(paramStaBinsType);

if strcmp(dateIntervalFlag, 'Month')
    dateList = datetimeList.Month;
elseif strcmp(dateIntervalFlag, 'Day')
    dateList = datetimeList;
end
dateTypes = unique(dateList);
dateTypesN = length(dateTypes);

paramStaBinsCount = zeros(paramStaBinsTypeN, 1) * nan;
dateIntervalCount = zeros(dateTypesN, 1) * nan;

lstArrayParamSta = zeros(paramStaBinsTypeN, dateTypesN) * nan;
for i = 1 : paramStaBinsTypeN
    paramStaBinIndex = (paramStaBins == paramStaBinsType(i));
    paramStaBinsCount(i) = sum(paramStaBinIndex);
    for j = 1 : dateTypesN
        dateListIndex =  (dateList == dateTypes(j));
        dateIntervalCount(j) = sum(dateListIndex);
        lstArrayParamSta(i, j) = mean(lstArray(paramStaBinIndex, dateListIndex), 'all', 'omitnan');
    end
end
paramStaBinsRange = [paramStaIntervals(paramStaBinsType)',  paramStaIntervals(paramStaBinsType+1)'];
end


% 冰川表面温度时间序列图。
function f = gstFigure(dateList, lstArrayParamSta, paramStabins, paramBinsCount, dateType, ...
    daynight, dataType, region, yearStr, topoType)
% 输入参数：
% dateList 是时间序列横轴。
% lstArrayParamSta 是冰川表面温度分级矩阵。
% bins 是分级范围。
% daynight 是昼夜标记。
% region 是区域。
% yearStr 是年份字符串。
% topoType 是用于分级的地形参数。
% 输出参数：
% f 是Figure对象。
f = figure; f.Position = [500 300 1500 500]; f.Visible = 'off';
if ismember(dataType, {'MOD11A1', 'MYD11A1'})
    lstArrayParamSta = lstArrayParamSta-273.15;
end
plot(dateList, lstArrayParamSta, 'o-',...
    dateList, zeros(size(lstArrayParamSta, 2), 1), '--k');
if strcmp(dateType, 'Month')
    dateType2 = 'monthly';
elseif strcmp(dateType, 'Day')
    dateType2 = 'daily';
end

if ismember(dataType, {'MOD11A1', 'MYD11A1'})
    title(sprintf('Time series of %s %stime glacier surface temperature (%s) in %s, %s (%s ranges)',...
        dateType2, daynight, dataType, region, yearStr, topoType));
    ylabel('Glacier surface temperature (°C)');
elseif ismember(dataType, {'MOD10A1', 'MYD10A1'})
    title(sprintf('Time series of %s glacier surface albedo (%s) in %s, %s (%s ranges)',...
        dateType2, dataType, region, yearStr, topoType));
    ylabel('Glacier surface Albedo');
end

xlabel(dateType);
xlim([1 12]);
% if strcmp(daynight, 'Day')
%     ylim([-30 40]);
% elseif strcmp(daynight, 'Night')
%     ylim([-45 5]);
% end

paramStabinsN = size(paramStabins, 1); lgdString = cell(paramStabinsN, 1);
if strcmp(topoType, 'Elevation')
    topoUnit = 'm';
elseif strcmp(topoType, 'Slope')
    topoUnit = '°';
end
for j = 1 : paramStabinsN
    lgdString{j} = sprintf('%s range: %d-%d %s, Count %d', topoType, paramStabins(j, 1),...
        paramStabins(j, 2), topoUnit, paramBinsCount(j));
end
legend(lgdString, 'Location', 'northeastoutside');
end
