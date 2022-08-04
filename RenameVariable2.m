
rootPath = 'G:\HMA_LST_Albedo\Data';
stepsPath = fullfile(rootPath, 'Basic\GlacierAreaInPixel');
modisArrayStepPath = fullfile(stepsPath, 'Step10_HMA_Array_Matlab');
modisMatList = dir(modisArrayStepPath);
modisMatList = {modisMatList(3:end).name}';
modisMatListN = length(modisMatList);
for i = 1 : modisMatListN
    modisMatPath = fullfile(modisArrayStepPath, modisMatList{i});
    variableList = who('-file', modisMatPath);
    if ismember({'dateList', 'paraMatrix'}, variableList)
        load(modisMatPath)
        modisDateList = dateList;
        modisMatrix = single(paraMatrix);
        save(modisMatPath, 'modisDateList', 'modisMatrix');
        disp(modisMatList{i});
    end
end
