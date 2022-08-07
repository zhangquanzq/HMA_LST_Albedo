%% HMA冰川的MODIS地表温度时间序列.

% !!! 坡度, 高程分级标准目前是等间距, 导致每个间隔内的像元样本数量差异大，看是否能改成保证每个间隔内数量一致
%   的分级标准. !!!

%% 标记和预设参数.
% 指定数据类型的标记. 1表示MOD10A1, 2表示MYD10A1, 3表示MOD11A1, 4表示MYD11A1.
flg1 = 1;
% 指定昼夜的标记. 1表示白天, 2表示晚上.
flg2 = 2;
% 指定面积比例的标记. 1表示100, 2表示95, 3表示90, 4表示85, 5表示80.
flg3 = 1;
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

seasons = struct('Winter0', 12, 'Winter1', [1, 2], 'Spring', [3, 4, 5], 'Summer', [6, 7, 8], ...
    'Autumn', [9, 10, 11]);

pctList = {'100', '95', '90', '85', '80'};
pct = pctList{flg3}; minPct = pctList{end};

topoList = {'Elev', 'Slp'};
unitList = {'meter', 'degree'};

% 高程和坡度区间.
elevEdges = 2000: 200: 8800;
slpEdges = 0: 10: 90;
topoEdges = {elevEdges, slpEdges};

% LST的有效qc取值.
qcValid = [0, 17];

%% 路径.
% 根目录.
rootDir = 'E:\HMA_LST_Albedo\Data';
stepsDir = fullfile(rootDir, 'GlacierAreaInPixel');
modisDir = fullfile(rootDir, dataType);

% 各输入数据存放的文件夹路径.
hmaPixelStaPctDir = fullfile(stepsDir, 'Step8_HMA_Pixel_Sta_Percent');
hmaPixelPctRasterDir = fullfile(stepsDir, 'Step9_HMA_Pixel_Percent_Raster');
modisMosaicDir = fullfile(modisDir, sprintf('%s_2_Mosaic%s_TIF', dataType, extent));

% 创建输出Mat文件的文件夹路径.
hmaMatDir = fullfile(stepsDir, 'Step10_HMA_Matlab');
if ~exist(hmaMatDir, 'dir')
    mkdir(hmaMatDir)
end

% 创建输出TIF文件的文件夹路径.
modisMeanDir = fullfile(modisDir, sprintf('%s_3_HMAGlacier%s_Mean_TIF', dataType, pct));
if ~exist(modisMeanDir, 'dir')
    mkdir(modisMeanDir);
end

% 创建存储各时间序列表格, 图的文件夹.
dataPctFolder = [dataType, '_HMAGlacier', pct];

tableDir = fullfile(stepsDir, 'Step11_HMA_Tables');
if ~exist(tableDir, 'dir')
    mkdir(tableDir)
end
tablePctDir = fullfile(tableDir, dataPctFolder);
if ~exist(tablePctDir, 'dir')
    mkdir(tablePctDir);
end

figDir = fullfile(stepsDir, 'Step12_HMA_Figures');
if ~exist(figDir, 'dir')
    mkdir(figDir)
end
figPctDir = fullfile(figDir, dataPctFolder);
if ~exist(figPctDir, 'dir')
    mkdir(figPctDir);
end
figAnnualDir = fullfile(figDir, 'Annual');
if ~exist(figAnnualDir, 'dir')
    mkdir(figAnnualDir)
end
figSeasonalDir = fullfile(figDir, 'Seasonal');
if ~exist(figSeasonalDir, 'dir')
    mkdir(figSeasonalDir)
end

% 最低冰川面积比例(80%)以上的像元FID编号的栅格文件.
hmaMinPctRasterName = ['HMA_ModisPixel_', cellsize, '_rgi60_FID_', minPct, 'percent.tif'];
hmaMinPctRasterPath = fullfile(hmaPixelPctRasterDir, hmaMinPctRasterName);

% 每个像元内高程, 坡度统计的矢量文件.
hmaMinPctStaElevName = ['HMA_ModisPixel_', cellsize, '_sta_elev_', minPct, 'percent.shp'];
hmaMinPctStaElevPath = fullfile(hmaPixelStaPctDir, hmaMinPctStaElevName);

hmaMinPctStaSlpName = ['HMA_ModisPixel_', cellsize, '_sta_slp_', minPct, 'percent.shp'];
hmaMinPctStaSlpPath = fullfile(hmaPixelStaPctDir, hmaMinPctStaSlpName);

% 每年的MODIS冰川像元值Mat文件名(用在sprintf函数中).
modisMatName = ['HMA_', dataType, '_%s_', minPct, 'percent.mat'];

% 按高程, 坡度分级统计的LST/Albedo时间序列Excel文件名(用在sprintf函数中).
lstXlsStr = [daynight, 'time_GST_timeseries_', dataType, '_', dateType, '_%s_%s.xlsx'];
albedoXlsStr = ['Albedo_timeseries_', dataType, '_', dateType, '_%s_%s.xlsx'];
lstFigStr = [daynight, 'time_GST_timeseries_', dataType, '_', dateType, '_%s_%s_%s.png'];
albedoFigStr = ['Albedo_timeseries_', dataType, '_', dateType, '_%s_%s_%s.png'];

% 获取最低冰川面积比例(80%)以上像元的位置索引和FID值.
[hmaMinPctLayer, hmaMinPctRef] = readgeoraster(hmaMinPctRasterPath);
minPctNodataValue = georasterinfo(hmaMinPctRasterPath).MissingDataIndicator;
geoTag = geotiffinfo(hmaMinPctRasterPath).GeoTIFFTags.GeoKeyDirectoryTag;

minPctIndexLayer = (hmaMinPctLayer ~= minPctNodataValue);
minPctPixelN = sum(minPctIndexLayer(:));
minPctFidList = hmaMinPctLayer(minPctIndexLayer);
[hmaRowN, hmaColN] = size(hmaMinPctLayer);

% 获取冰川MODIS像元的高程, 坡度统计表格, 按MODIS Pixel FID对属性表排序.
[~, staElevTable] = shaperead(hmaMinPctStaElevPath);
staElevTable = staElevTable(minPctFidList + 1);  % FID从0开始，加1从1开始.
elevMeanRecords = [staElevTable.MEAN]';

[~, staSlpTable] = shaperead(hmaMinPctStaSlpPath);
staSlpTable = staSlpTable(minPctFidList + 1);  % FID从0开始，加1从1开始.
slpMeanRecords = [staSlpTable.MEAN]';

topoMeanRecords = {elevMeanRecords, slpMeanRecords};

% 获取满足冰川面积百分比的像元编号.
if strcmp(pct, '100'); pctValue = 0.9999996427; else; pctValue = str2double(pct) / 100; end
validPctIndex = [staElevTable.Area_Pct]' >= pctValue;

% 获取冰川分区名称列表.
regionRecords = {staElevTable.FULL_NAME}';
regionList = unique(regionRecords);
regionListN = length(regionList);

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
    modisYearDir = fullfile(modisMosaicDir, modisYearFolder);
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
    break
    yearStr = num2str(yearList(i));
    
    % 从Mat文件中读取MODIS LST或Albedo数据, 并做质量控制.
    fmtStr = {yearStr, [daynight, '_', yearStr]};
    modisMatPath = fullfile(hmaMatDir, sprintf(modisMatName, fmtStr{round(flg1/2)}));
    [modisMatrix, modisDateList] = loadModisMatrix(modisMatPath, dataName, validPctIndex, qcValid);

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
        modisYearMeanLayer = zeros(hmaRowN, hmaColN, 'single') * nan;
        modisYearMeanLayer(minPctIndexLayer) = mean(modisMatrix, 2, 'omitnan');
        geotiffwrite(modisYearMeanPath, modisYearMeanLayer, hmaMinPctRef, ...
            GeoKeyDirectoryTag=geoTag, tifftags=struct('Compression','LZW'));
        disp(['输出: ', modisYearMeanName]);
    end
    
    % MODIS冰川区LST/Albedo的月平均值.
    monthList = modisDateList.Month;
    monthType = unique(monthList);
    modisMonthMeanLayer = zeros(hmaRowN, hmaColN, 'single') * nan;
    for j = 1 : length(monthType)
        monthNum = monthType(j); monthStr = num2str(monthNum, '%02d');
        modisMonthMeanName = replace(modisYearMeanName, yearStr, [yearStr, monthStr]);
        modisMonthMeanPath = fullfile(modisYearDir, modisMonthMeanName);
        if exist(modisMonthMeanPath, 'file')
            continue
        end
        modisMonthMeanVector = mean(modisMatrix(:, monthList == monthNum), 2, 'omitnan');
        modisMonthMeanLayer(minPctIndexLayer) = modisMonthMeanVector;
        geotiffwrite(modisMonthMeanPath, modisMonthMeanLayer, hmaMinPctRef, ...
            GeoKeyDirectoryTag=geoTag, tifftags=struct('Compression','LZW'));
        disp(['输出: ', modisMonthMeanName]);
    end
end

%% 按年周期创建不同分区冰川表面温度的时间序列图, 按高程, 坡度统计.
hmaYearMeanVector = zeros(yearListN, 1);
for i = 1 : yearListN
    break
    yearStr = num2str(yearList(i));
    
    % 从Mat文件中读取MODIS LST或Albedo数据, 并做质量控制.
    fmtStr = {yearStr, [daynight, '_', yearStr]};
    modisMatPath = fullfile(hmaMatDir, sprintf(modisMatName, fmtStr{round(flg1/2)}));
    [modisMatrix, modisDateList] = loadModisMatrix(modisMatPath, dataName, validPctIndex, qcValid);

    % HMA地区整体年均值.
    hmaYearMeanVector(i) = mean(modisMatrix, 'all', 'omitnan');

    % 分区统计.
    for j =  1 : regionListN
        regionName1 = regionList{j};
        regionName2 = replace(regionName1, '/', '_');

        % 创建存储各时间序列图的文件夹.
        figRegionDir = fullfile(figPctDir, regionName2);
        if ~exist(figRegionDir, 'dir')
            mkdir(figRegionDir);
        end

        % 将按高程, 坡度分级统计的结果存储为Excel文件, 并制图.
        for k = 1: length(topoList)
            topo = topoList{k}; edges = topoEdges{k}; unit = unitList{k};

            % 分级.
            regionIndex = strcmp(regionName1, regionRecords);
            subModisMatrix = modisMatrix(regionIndex, :);
            subTopoMeanRecords = topoMeanRecords{k}(regionIndex);
            [topoBinsRange, topoBinsCount, topoDateTypes, topoDateCount, topoMeanMatrix] = ...
                intervalMean(subModisMatrix, subTopoMeanRecords, edges, modisDateList, ...
                dateType);

            % 存Excel.
            lstXlsName = sprintf(lstXlsStr, regionName2, topo);
            albedoXlsName = sprintf(albedoXlsStr, regionName2, topo);
            xlsName = {albedoXlsName, lstXlsName};
            xlsPath = fullfile(tablePctDir, xlsName{round(flg1/2)});
            xlsExist = exist(xlsPath, 'file');
            if ~xlsExist || (xlsExist && ~ismember(yearStr, sheetnames(xlsPath)))
                topoBinsRangeN = length(topoBinsCount); topoDateTypeN = length(topoDateTypes);
                topoBinsRangeCell = cell(topoBinsRangeN, 1);
                for m = 1 : topoBinsRangeN
                    topoBinsMin = topoBinsRange(m, 1); topoBinsMax = topoBinsRange(m, 2);
                    topoBinsRangeCell(m) = {sprintf('%d-%d %s', topoBinsMin, topoBinsMax, unit)};
                end
                topoStaCell = cell(topoBinsRangeN + 2, topoDateTypeN + 2);
                topoStaCell(1, 2) = {dateType};
                topoStaCell(1, 3:end) = num2cell(topoDateTypes');
                topoStaCell(2, 1) = {sprintf('%s range', topo)};
                topoStaCell(2, 2) = {'Count'};
                topoStaCell(2, 3:end) = num2cell(topoDateCount');
                topoStaCell(3:end, 1) = topoBinsRangeCell;
                topoStaCell(3:end, 2) = num2cell(topoBinsCount);
                topoStaCell(3:end, 3:end) = num2cell(topoMeanMatrix);
                writecell(topoStaCell, xlsPath, 'Sheet', yearStr)
                fprintf('添加%s年统计值到%s\n', yearStr, xlsName{round(flg1/2)})
            end

            % 制图.
            LstFigName = sprintf(lstFigStr, regionName2, yearStr, topo);
            albedoFigName = sprintf(albedoFigStr, regionName2, yearStr, topo);
            figName = {albedoFigName, LstFigName};
            figPath = fullfile(figRegionDir, figName{round(flg1/2)});
            if ~exist(figPath, 'file')
                fig = tsFigure(topoDateTypes, topoMeanMatrix, topoBinsRange, topoBinsCount, ...
                    dateType, daynight, dataType, regionName1, yearStr, topo);
                %  exportgraphics(fig, figPath)
                print(fig, figPath, '-dpng', '-r200')
                fprintf('输出分级时间序列图 %s\n', figPath)
                close all
            end
        end        
    end
end

% 年纪变化趋势图.
fig = figure; fig.Position = [400 300 1000 400];
if strcmp(dataName, 'LST')
    hmaYearMeanVector = hmaYearMeanVector - 273.15;
end

p2 = polyfit(yearList, hmaYearMeanVector, 1);
y1 = p2(1) * 2000 + p2(2); y2 = p2(1) * 2020 + p2(2);
plot(yearList, hmaYearMeanVector, 'o-', [2000, 2020], [y1, y2], MarkerFaceColor='b');
% legend(dataName, 'Trend', Location='bestoutside', Orientation='horizontal')
xlabel('Year');

if strcmp(dataName, 'LST')
    ylabel('Glacier surface temperature (°C)')

    titleStr = ['Annual Mean of %s-time GST in HMA from %d to %d\n' ...
    '(%s, Glacier pixel threshold: %s%%)'];
    t = title(sprintf(titleStr, daynight, yearList(1), yearList(end), dataType, pct));
    t.FontWeight = 'bold';

    figNameStr = 'Annual Mean of %s-time GST in HMA from %d to %d, %s, %s%%.png';
    figName = sprintf(figNameStr, daynight, yearList(1), yearList(end), dataType, pct);
elseif strcmp(dataName, 'Albedo')
    ylabel('Glacier surface Albedo')

    titleStr = ['Annual Mean of Glacier surface Albedo in HMA from %d to %d\n' ...
        '(%s, Glacier pixel threshold: %s%%)'];
    t = title(sprintf(titleStr, yearList(1), yearList(end), dataType, pct));
    t.FontWeight = 'bold';

    figNameStr = 'Annual Mean of Glacier surface Albedo in HMA from %d to %d, %s, %s%%.png';
    figName = sprintf(figNameStr, yearList(1), yearList(end), dataType, pct);
end

figPath = fullfile(figAnnualDir, figName);
if ~exist(figPath, 'file')
%     exportgraphics(fig, figPath)
    print(fig, figPath, '-dpng', '-r200')
    fprintf('输出年际时间序列图 %s\n', figPath)
end
close all

%% 季节和年纪变化趋势图.
% 季节变化趋势数据. 冬季为上一年度的12月和本年度的1, 2月.
modisSeasonMeanMatrix = zeros(4, yearListN);
for i = 2 : yearListN
    % 上一年.
    if i ~= 1
        yearStr0 = num2str(yearList(i-1)); fmtStr0 = {yearStr0, [daynight, '_', yearStr0]};
        matPath0 = fullfile(hmaMatDir, sprintf(modisMatName, fmtStr0{round(flg1/2)}));
        [modisMatrix0,modisDateList0] = loadModisMatrix(matPath0,dataName,validPctIndex,qcValid);
        modisMonthList0 = modisDateList0.Month;
    end

    % 当前年.
    yearStr1 = num2str(yearList(i)); fmtStr1 = {yearStr1, [daynight, '_', yearStr1]};
    matPath1 = fullfile(hmaMatDir, sprintf(modisMatName, fmtStr1{round(flg1/2)}));
    [modisMatrix1,modisDateList1] = loadModisMatrix(matPath1,dataName,validPctIndex,qcValid);
    modisMonthList1 = modisDateList1.Month;

    % 若是第一年, 缺上一年12月数据, 冬季仅包括当前年1, 2月, 否则冬季包括上一年12月和当前年1, 2月.
    if i == 1
        winterMatrix = modisMatrix1(:, ismember(modisMonthList1, seasons.Winter1));
    else
        winter0Matrix = modisMatrix0(:, ismember(modisMonthList0, seasons.Winter0));
        winter1Matrix = modisMatrix1(:, ismember(modisMonthList1, seasons.Winter1));
        winterMatrix = [winter1Matrix, winter0Matrix];
    end
    springMatrix = modisMatrix1(:, ismember(modisMonthList1, seasons.Spring));
    summerMatrix = modisMatrix1(:, ismember(modisMonthList1, seasons.Summer));
    autumnMatrix = modisMatrix1(:, ismember(modisMonthList1, seasons.Autumn));
    
    winterMean = mean(winterMatrix, 'all', 'omitnan');
    sprintMean = mean(springMatrix, 'all', 'omitnan');
    summerMean = mean(summerMatrix, 'all', 'omitnan');
    autumnMean = mean(autumnMatrix, 'all', 'omitnan');
    modisSeasonMeanMatrix(:, i) = [winterMean, sprintMean, summerMean, autumnMean];
end

% 季节变化趋势图.
fig = figure; fig.Position = [400 300 1000 400];
if strcmp(dataName, 'LST')
    modisSeasonMeanMatrix = modisSeasonMeanMatrix - 273.15;
end
for i = 1 : 4
    modisSeasonMeanVector = modisSeasonMeanMatrix(i, :);
    p1 = polyfit(yearList, modisSeasonMeanVector, 1);
    y1 = p1(1) * 2000 + p1(2); y2 = p1(1) * 2020 + p1(2);
    plot(yearList, modisSeasonMeanVector, 'o-', [2000, 2020], [y1, y2], 'k')
    hold on
end
legend('Winter', '', 'Spring', '', 'Summer', '', 'Autumn', '', Location='bestoutside', ...
    Orientation='horizontal')
xlabel('Year')

if strcmp(dataName, 'LST')
    ylabel('Glacier surface temperature (°C)');

    titleStr = ['Seasonal Mean of %s-time Glacier surface temperature in HMA from %d to %d\n' ...
        '(%s, Glacier pixel threshold: %s%%)'];
    t = title(sprintf(titleStr, daynight, yearList(1), yearList(end), dataType, pct));
    t.FontWeight = 'bold';

    figNameStr = 'Seasonal Mean of %s-time GST in HMA from %d to %d, %s, %s%%.png';
    figName = sprintf(figNameStr, daynight, yearList(1), yearList(end), dataType, pct);

elseif strcmp(dataName, 'Albedo')
    ylabel('Glacier surface Albedo');
    titleStr = ['Seasonal Mean of Glacier surface Albedo in HMA from %d to %d\n' ...
        '(%s, Glacier pixel threshold: %s%%)'];
    t = title(sprintf(titleStr, yearList(1), yearList(end), dataType, pct));
    t.FontWeight = 'bold';

    figNameStr = 'Seasonal Mean of Glacier surface Albedo in HMA from %d to %d, %s, %s%%.png';
    figName = sprintf(figNameStr, yearList(1), yearList(end), dataType, pct);
end

figPath = fullfile(figSeasonalDir, figName);
if ~exist(figPath, 'file')
%     exportgraphics(fig, figPath)
    print(fig, figPath, '-dpng', '-r200')
    fprintf('输出季节时间序列图 %s\n', figPath)
end
close all

%% 自定义函数.
% 从Mat文件中读取数据质量控制后的MODIS LST/Albedo数据.
function [modisMatrix, modisDateList] = loadModisMatrix(matPath, dataName, validPctIndex, qcValid)
load(matPath, 'modisMatrix', 'modisDateList');
if strcmp(dataName, 'LST')
    load(matPath, 'qcMatrix');
    modisMatrix(~ismember(qcMatrix, qcValid)) = nan;
elseif strcmp(dataName, 'Albedo')
    modisMatrix(modisMatrix == 0 | modisMatrix > 100) = nan;
    modisMatrix = modisMatrix ./ 100;
end
modisMatrix(~validPctIndex, :) = nan;  % 保留指定冰川面积比例的像元.
end

% 按某一参数的统计值划分区间, 计算每一区间的均值.
function [binsRange, binsCount, dateTypes,  dateCount, meanMatrix] = ...
    intervalMean(dataMatrix, staRecords, edges, datetimeList, dateType)
% 输入参数:
% dataMatrix 是冰川表面温度, 反照率值矩阵, 行表示不同的空间位置, 列表示日期(年内).
% staRecords 是参数统计值列向量, 每个值表示dataMatrix中每行对应的像元位置的参数统计.
% edges 是参数统计值的划分区间. 例: 坡度取值0-90度, 可划分为5度间隔的多个区间(0: 5: 90).
% datetimeList 是某年内的日期列表.
% dateType 是日期划分区间标记。可取'Month'和'Day'.

% 输出参数:
% binsRange 是参数统计值区间范围矩阵, 每行的2数字表示区间值域, 行数为参数统计值区间的个数.
% binsCount 是参数统计值区间内的统计值个数向量, 行数为参数统计值区间的个数.
% dateTypes 是日期列表按dateType划分后的日期区间值.
% dateCount 是每个日期区间包含的日期个数.
% meanMatrix 是按参数取值区间计算的地表温度值矩阵.

% 注意: staRecords中相邻两个统计值的所属区间没有特定顺序.

staBins = discretize(staRecords, edges);
staBinsType = unique(staBins);
staBinsTypeN = length(staBinsType);

if strcmp(dateType, 'Month')
    dateList = datetimeList.Month;
elseif strcmp(dateType, 'Day')
    dateList = datetimeList;
end
dateTypes = unique(dateList);
dateTypesN = length(dateTypes);

binsCount = zeros(staBinsTypeN, 1) * nan;
dateCount = zeros(dateTypesN, 1) * nan;
meanMatrix = zeros(staBinsTypeN, dateTypesN) * nan;
for i = 1 : staBinsTypeN
    staBinIndex = (staBins == staBinsType(i));
    binsCount(i) = sum(staBinIndex);
    for j = 1 : dateTypesN
        dateListIndex = (dateList == dateTypes(j));
        dateCount(j) = sum(dateListIndex);
        meanMatrix(i, j) = mean(dataMatrix(staBinIndex, dateListIndex), 'all', 'omitnan');
    end
end
binsRange = [edges(staBinsType)', edges(staBinsType+1)'];
end

% 冰川表面温度, 反照率时间序列图.
function f = tsFigure(dateList, meanMatrix, binsRange, binsCount, dateType, daynight, dataType, ...
    region, yearStr, topoType)
% 输入参数:
% dateList 时间序列横轴.
% meanMatrix 冰川表面温度分级矩阵.
% binsRange 分级范围.
% binsCount 每级别的记录数量.
% dateType 日期统计间隔.
% daynight 昼夜标记.
% dataType 数据类型.
% region 冰川分区.
% yearStr 年份字符串.
% topoType 用于分级的地形参数.

% 输出参数:
% f Figure对象.

f = figure; f.Position = [500 300 1500 500]; f.Visible = 'off';
if ismember(dataType, {'MOD11A1', 'MYD11A1'})
    meanMatrix = meanMatrix - 273.15;
end
plot(dateList, meanMatrix, 'o-', ...
    dateList, zeros(size(meanMatrix, 2), 1), '--k');
if strcmp(dateType, 'Month')
    dateType2 = 'monthly';
elseif strcmp(dateType, 'Day')
    dateType2 = 'daily';
end

if ismember(dataType, {'MOD11A1', 'MYD11A1'})
    title(sprintf('Timeseries of %s %stime GST (%s) in %s, %s (%s ranges)', dateType2, daynight, ...
        dataType, region, yearStr, topoType))
    ylabel('Glacier surface temperature (°C)')
elseif ismember(dataType, {'MOD10A1', 'MYD10A1'})
    title(sprintf('Timeseries of %s glacier surface albedo (%s) in %s, %s (%s ranges)', ...
        dateType2, dataType, region, yearStr, topoType))
    ylabel('Glacier surface Albedo')
end

xlabel(dateType); xlim([1 12]);
% if strcmp(daynight, 'Day')
%     ylim([-30 40]);
% elseif strcmp(daynight, 'Night')
%     ylim([-45 5]);
% end

stabinsN = size(binsRange, 1); lgdString = cell(stabinsN, 1);
if strcmp(topoType, 'Elev')
    topoUnit = 'm';
elseif strcmp(topoType, 'Slp')
    topoUnit = '°';
end
for j = 1 : stabinsN
    lgdString{j} = sprintf('%s range: %d-%d %s, Count %d', topoType, binsRange(j, 1),...
        binsRange(j, 2), topoUnit, binsCount(j));
end
legend(lgdString, Location='northeastoutside')
end
