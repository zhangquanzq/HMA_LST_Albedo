%% 标记和预设参数。
% 指定数据类型的标记。1表示MOD10A1，2表示MYD10A1，3表示MOD11A1，4表示MYD11A1。
flg1 = 1;
% 指定昼夜的标记。1表示白天，2表示晚上。
flg2 = 1;
% 指定面积比例的标记。1表示100，2表示95，3表示90，4表示85，5表示80。
flg3 = 1;

dataTypeList = {'MOD10A1', 'MYD10A1', 'MOD11A1', 'MYD11A1'};
daynightList = {'Day', 'Night'};
percentList = {'100', '95', '90', '85', '80'};
yearListList = {2002 : 2019, 2000 : 2019};


dataType = dataTypeList{flg1};
daynight = daynightList{flg2};
percent = percentList{flg3};
yearList = yearListList{mod(flg1, 2) + 1};
yearListN = length(yearList);

rootPath = 'F:\HMA_LST_Albedo\Data';
stepsPath = fullfile(rootPath, 'Basic\GlacierAreaInPixel');
modisPath = fullfile(rootPath, dataType);

modisArrayStepPath = fullfile(stepsPath, 'Step10_HMA_Array_Matlab - 副本');
dataPercentTifFolder = [dataType, '_3_HMAGlacier', percent, '_TIF'];
for i = 1 : yearListN
    yearStr = num2str(yearList(i));
    dataYearFolder = [dataType, '_', yearStr, 'XXX_TIF'];
    modisYearFolderPath = fullfile(modisPath, dataPercentTifFolder, dataYearFolder);

    % 读取MODIS LST/Albedo数据的文件名列表。
    paraYearDateTifName = {[dataType, '*Albedo.tif'], [dataType, '*LST_', daynight, '.tif']};
    paraYearDateTifName = paraYearDateTifName{round(flg1/2)};
    paraYearDateTifList = dir(fullfile(modisYearFolderPath, paraYearDateTifName));
    paraYearDateTifList = {paraYearDateTifList.name}';
    paraYearDateTifListN = length(paraYearDateTifList);

    % 将每一年HMA冰川区MODIS数据存储为Mat文件，提高程序运行效率。
    paraMatName = {sprintf('hma_%s_%s_%spercent.mat', dataType, yearStr, percent), ...
        sprintf('hma_%s_%s_%s_%spercent.mat', dataType, daynight, yearStr, percent)};
    paraMatName = paraMatName{round(flg1/2)};
    paraMatPath = fullfile(modisArrayStepPath, paraMatName);
    load(paraMatPath, 'paraMatrix')
    delete(paraMatPath)

    % 从MODIS LST/Albedo数据的文件名中获取数据日期列表。
    dateList = cell(paraYearDateTifListN, 1);
    for j = 1 : paraYearDateTifListN
        lstTifName = paraYearDateTifList{j};
        dateList{j} = yday2ymd(lstTifName(10:16));
    end
    dateList = datetime(dateList, 'InputFormat', 'yyyy/MM/dd');

    save(paraMatPath, 'paraMatrix', 'dateList');
    disp(matfile(paraMatPath))
end