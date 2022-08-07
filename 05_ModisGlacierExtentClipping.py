# coding=utf-8
# 按照不同的冰川占像元比例裁剪HMA地区的MODIS LST/Albedo数据。
import arcpy
import os

# 指定数据类型的标记。1表示MOD10A1，2表示MYD10A1，3表示MOD11A1，4表示MYD11A1。
flg1 = 4
# 确定裁剪的MODIS像元中，冰川面积比例。 1表示100%，2表示95%，3表示90%，4表示85%, 5表示80%。
flg2 = 5

modisType = ['MOD10A1', 'MYD10A1', 'MOD11A1', 'MYD11A1'][flg1 - 1]
region = ['HMA', 'HMA', 'CN', 'CN'][flg1 - 1]
resolution = ['500m', '500m', '1km', '1km'][flg1 - 1]
glacierAreaPercent = [100, 95, 90, 85, 80][flg2 - 1]

arcpy.env.overwriteOutput = True
# 输入输出路径。
rootPath = 'H:\\HMA_LST_Albedo\\Data\\'

modisMosaicFolder = '{0}_2_Mosaic{1}_TIF'.format(modisType, region)
modisMosaicPath = os.path.join(rootPath, modisType, modisMosaicFolder)

modisGlacierFolder = '{0}_3_HMAGlacier{1}_TIF'.format(modisType, glacierAreaPercent)
modisGlacierPath = os.path.join(rootPath, modisType, modisGlacierFolder)
if not os.path.exists(modisGlacierPath):
    os.makedirs(modisGlacierPath)

hmaModisPixelTifName = 'HMA_ModisPixel_{0}_rgi60_{1}percent.tif'.format(resolution, glacierAreaPercent)
hmaModisPixelTifFolderPath = os.path.join(rootPath, 'Basic\\GlacierAreaInPixel\\Step9_HMA_Pixel_Percent_Raster')
hmaModisPixelTifPath = os.path.join(hmaModisPixelTifFolderPath, hmaModisPixelTifName)

# 裁剪MODIS LST影像中的冰川覆盖像元。
arcpy.env.workspace = modisMosaicPath
modisYearPathList = arcpy.ListWorkspaces('*', 'Folder')
modisErrorFileList = []
for modisYearPath in modisYearPathList:
    modisYearFolder = os.path.basename(modisYearPath)
    modisGlacierYearPath = os.path.join(modisGlacierPath, modisYearFolder)
    if not os.path.exists(modisGlacierYearPath):
        os.makedirs(modisGlacierYearPath)
    arcpy.env.workspace = modisYearPath
    modisTifList = arcpy.ListRasters('*', 'TIF')
    for modisTif in modisTifList:
        modisGlacierTifPath = os.path.join(modisGlacierYearPath, modisTif)
        if not os.path.exists(modisGlacierTifPath):  # 如果文件不存在，裁剪生成特定比例的MODIS冰川像元。
            if resolution == '500m':
                modisQcRaster = arcpy.sa.SetNull(modisTif, modisTif, 'VALUE > 100')
            else:
                modisQcRaster = arcpy.Raster(modisTif)
            outTimes = arcpy.sa.Times(modisQcRaster, hmaModisPixelTifPath)
            outTimes.save(modisGlacierTifPath)
            print(u'{0} 完成HMA地区冰川面积占比为 {1}% 的MODIS像元提取.'.format(modisTif, glacierAreaPercent))
        else:  # 如果文件存在，判断文件是否正常。将不正常的文件路径单独存储，以备后续检查。
            try:
                modisRaster = arcpy.Raster(modisGlacierTifPath)
            except:
                modisErrorFileList.append(modisGlacierTifPath + '\n')
                print(u'{0} 文件错误.'.format(modisTif))
errorListFile = os.path.join(modisMosaicPath, '{0} {1}% error files.txt'.format(modisType, glacierAreaPercent))
with open(errorListFile, 'w') as f:
    f.writelines(modisErrorFileList)
