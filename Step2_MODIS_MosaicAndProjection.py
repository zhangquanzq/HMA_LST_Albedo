# %% coding=utf-8
# 将分年份, 日期组织的MODIS地表温度和反照率HDF文件拼接成一幅TIF文件, 然后投影转换.
# 注意: 此程序拼接得到的TIF文件中, 原HDF为nodata的数值变成数字, 需要下一步将其转为Nodata.
import arcpy
import os
import sys
from collections import Counter
from glob import glob
from os import path

'''
MxD10A1数据层: ['Snow_cover','Basic_QA','Algorithm_QA','NDSI','Albedo','Orbit_pnt','Granule_pnt']
MxD11A1数据层: ['LST_Day','QC_Day','Time_Day','Angle_Day','LST_Night','QC_Night','Time_Night',
               'Angle_Night', 'Emis_31', 'Emis_32', 'Clear_Day', 'Clear_Night']
'''

# 指定MODIS数据类型的标记. 1表示MOD10A1, 2表示MYD10A1, 3表示MOD11A1, 4表示MYD11A1.
flg1 = 4
# 指定区域的标记. 1表示HMA, 2表示CN.
flg2 = 2

# LST和Albedo, 以及MOD和MYD的属性.
lstPixelType, albedoPixelType = '16_BIT_UNSIGNED', '8_BIT_UNSIGNED'
lstLayerList, albedoLayerList = ['LST_Day', 'QC_Day', 'LST_Night', 'QC_Night', 'Emis_31',
                                 'Emis_32'], ['Albedo']
lstLayerNumList, albedoLayerNumList = ['0', '1', '4', '5', '8', '9'], ['4']
modYearList, mydYearList = range(2000, 2020), range(2002, 2020)

# 四种MODIS数据的属性字典.
mod10a1Atts = {'modisType': 'MOD10A1', 'pixelType': albedoPixelType, 'yearList': modYearList,
               'layerList': albedoLayerList, 'layerNumList': albedoLayerNumList}
myd10a1Atts = {'modisType': 'MYD10A1', 'pixelType': albedoPixelType, 'yearList': mydYearList,
               'layerList': albedoLayerList, 'layerNumList': albedoLayerNumList}
mod11a1Atts = {'modisType': 'MOD11A1', 'pixelType': lstPixelType, 'yearList': modYearList,
               'layerList': lstLayerList, 'layerNumList': lstLayerNumList}
myd11a1Atts = {'modisType': 'MYD11A1', 'pixelType': lstPixelType, 'yearList': mydYearList,
               'layerList': lstLayerList, 'layerNumList': lstLayerNumList}

# 预设参数.
modisAtts = [mod10a1Atts, myd10a1Atts, mod11a1Atts, myd11a1Atts][flg1 - 1]
region = ['HMA', 'CN'][flg2 - 1]
modisType = modisAtts['modisType']
pixelType = modisAtts['pixelType']
yearList = modisAtts['yearList']
layerList = modisAtts['layerList']
layerNumList = modisAtts['layerNumList']
layerListN = len(layerList)

# Albedo数据只处理HMA区.
if modisType in ['MOD10A1', 'MYD10A1'] and region != 'HMA':
    print('Albedo数据只能处理HMA区. 当前数据为{0}, 区域为{1}, 程序退出.'.format(modisType, region))
    sys.exit()

# ArcPy环境参数设置.
arcpy.env.overwriteOutput = True
arcpy.env.pyramid = 'NONE'
# arcpy.env.rasterStatistics = 'NONE'
# arcpy.env.nodata = 'NONE'

# %% 路径.
# 根目录.
rootPath = r'G:\HMA_LST_Albedo'
dataPath = path.join(rootPath, 'Data', modisType)
# 输入数据路径.
modisTilesPath = path.join(dataPath, '{0}_1_Tile{1}_HDF'.format(modisType, region))
# 输出数据路径.
modisMosaicPath = path.join(dataPath, '{0}_2_Mosaic{1}_TIF'.format(modisType, region))
if not path.exists(modisMosaicPath):
    os.makedirs(modisMosaicPath)
modisPrjPath = path.join(dataPath, '{0}_2_Prj{1}_TIF'.format(modisType, region))
if not path.exists(modisPrjPath):
    os.makedirs(modisPrjPath)
# 空间坐标系统文件路径.
prjPath = path.join(dataPath, 'Modis_Sinusoidal.prj')

# %% 数据处理.
# 按年份拼接MODIS数据图层.
for year in yearList:
# for year in range(2021, 2022):
    modisTilesYearPath = path.join(modisTilesPath, '{0}_{1}XXX_HDF'.format(modisType, year))
    modisMosaicYearPath = path.join(modisMosaicPath, '{0}_{1}XXX_TIF'.format(modisType, year))
    if not path.exists(modisMosaicYearPath):
        os.makedirs(modisMosaicYearPath)

    # 获取已完成所有数据层拼接的日期列表.
    arcpy.env.workspace = modisMosaicYearPath
    modisMosaicLayerDateList, modisMosaicCompleteDateList = [], []
    for layerName in layerList:
        mosaicLayerList = arcpy.ListRasters('{0}*{1}.tif'.format(modisType, layerName))
        for mosaicLayer in mosaicLayerList:
            modisMosaicLayerDateList.append(mosaicLayer.split('.')[1][1:])
    modisMosaicDateCounter = Counter(modisMosaicLayerDateList)
    modisMosaicDateNumberList = modisMosaicDateCounter.values()
    modisMosaicDateStrList = modisMosaicDateCounter.keys()
    for i in range(len(modisMosaicDateNumberList)):
        if modisMosaicDateNumberList[i] == layerListN:
            modisMosaicCompleteDateList.append(modisMosaicDateStrList[i])

    # 提取MODIS数据层并拼接它们, 跳过已完成所有图层拼接的日期.
    for modisTilesYearDayPath in glob(path.join(modisTilesYearPath, '*')):
        yearDay = path.basename(modisTilesYearDayPath).split('_')[-1]
        if yearDay in modisMosaicCompleteDateList:
            continue  # 跳过已完成所有数据层拼接的日期.

        # 从MODIS HDF文件中提取数据层，并保存到临时文件夹中.
        print(u'提取{0} {1} {2}'.format(region, yearDay, modisType))
        tmpPath = path.join(rootPath, 'tmp_{0}_{1}'.format(modisType, yearDay))
        if not path.exists(tmpPath):
            os.mkdir(tmpPath)
        for modisTilePath in glob(path.join(modisTilesYearDayPath, '*.hdf')):
            for i in range(layerListN):
                layerName = layerList[i]
                modisMosaicFileName = '{0}.A{1}.{2}.tif'.format(modisType, yearDay, layerName)
                modisMosaicFilePath = path.join(modisMosaicYearPath, modisMosaicFileName)
                if not arcpy.Exists(modisMosaicFilePath):
                    modisTileName = path.basename(modisTilePath)
                    modisLayerName = modisTileName.replace('.hdf', '_{0}.tif'.format(layerName))
                    modisLayerPath = path.join(tmpPath, modisLayerName)
                    if not arcpy.Exists(modisLayerPath):
                        try:
                            arcpy.ExtractSubDataset_management(modisTilePath, modisLayerPath,
                                                               layerNumList[i])
                        except:
                            errorHDFText = path.join(rootPath, modisType,
                                                     'errorHDF{0}.txt'.format(yearDay))
                            with open(errorHDFText, 'a') as f:
                                f.write(modisTilePath + '\n')

        # 拼接提取的数据层.
        print(u'拼接{0} {1} {2}'.format(region, yearDay, modisType))
        for layerName in layerList:
            arcpy.env.workspace = tmpPath
            modisTileList = arcpy.ListRasters("{0}*_{1}.tif".format(modisType, layerName))
            modisMosaicFileName = '{0}.A{1}.{2}.tif'.format(modisType, yearDay, layerName)
            modisMosaicFilePath = path.join(modisMosaicYearPath, modisMosaicFileName)
            modisTileStr = ';'.join(modisTileList)
            if not arcpy.Exists(modisMosaicFilePath):
                gdbName = '{0}_{1}.gdb'.format(modisType, yearDay)
                gdbPath = path.join(rootPath, gdbName)
                mosaicPath = path.join(gdbPath, layerName)
                arcpy.CreateFileGDB_management(rootPath, gdbName)
                arcpy.CreateMosaicDataset_management(gdbPath, layerName, prjPath)
                arcpy.AddRastersToMosaicDataset_management(mosaicPath, 'Raster Dataset', tmpPath,
                                                           filter='*{0}.tif'.format(layerName))
                # 以下pixel_type和nodata_value仅是Emis的参数设置. , nodata_value='255'
                arcpy.CopyRaster_management(mosaicPath, modisMosaicFilePath, pixel_type=pixelType)
                if 'LST' in path.basename(modisMosaicFilePath):
                    arcpy.SetRasterProperties_management(modisMosaicFilePath, nodata='1 0')
                arcpy.Delete_management(gdbPath)

        # 删除临时文件夹.
        arcpy.Delete_management(tmpPath)

# 投影转换.
for year in yearList:
    modisMosaicYearPath = path.join(modisMosaicPath, '{0}_{1}XXX_TIF'.format(modisType, year))
    modisPrjYearPath = path.join(modisPrjPath, '{0}_{1}XXX_TIF'.format(modisType, year))
    if not path.exists(modisPrjYearPath):
        os.makedirs(modisPrjYearPath)
    modisMosaicTifPathList = glob(path.join(modisMosaicYearPath, '*.tif'))
    for modisMosaicFilePath in modisMosaicTifPathList:
        modisPrjTifPath = modisMosaicFilePath.replace('MosaicCN', 'PrjCN')
        if not path.exists(modisPrjTifPath):
            print(modisPrjTifPath)
            dirName = r'F:\AMSR_MODIS_Fusion\Data\AMSR2_2_CN_TIF\L3.TB7GHz_10\2014\03'
            snapRasterName = 'GW1AM2_20140300_01M_EQMA_L3SGT07HA2220220_BtH.tif'
            arcpy.env.snapRaster = path.join(dirName, snapRasterName)
            if '2014151' in modisMosaicFilePath:
                arcpy.ProjectRaster_management(modisMosaicFilePath, modisPrjTifPath,
                                               arcpy.SpatialReference(4326), 'NEAREST', '0.01')
