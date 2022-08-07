# coding=utf-8
# HMA冰川范围MODIS像元内的高程, 坡度统计.
import arcpy
from os import path

# 标记和预设参数.
# 指定MODIS数据类型的标记. 1表示MYD11A1, 2表示MOD11A1, 3表示MYD10A1, 4表示MOD10A1.
flg1 = 3
# 指定地形参数的标记. 1表示高程, 2表示坡度.
flg2 = 1

# 预设参数.
modisType = ['MYD11A1', 'MOD11A1', 'MYD10A1', 'MOD10A1'][flg1 - 1]
cellsize = ['1km', '1km', '500m', '500m'][flg1 - 1]
topoType = ['elev', 'slp'][flg2 - 1]
staTypeList = ['Maximum', 'Mean', 'Minimum', 'Range', 'Std']

# 路径.
rootDir = r'G:\HMA_LST_Albedo'
dataDir = path.join(rootDir, 'Data')

srtmPath = path.join(dataDir, r'Basic\Raster', 'HMA_SRTM_3arc_' + topoType + '.tif')

stepsDir = path.join(dataDir, 'GlacierAreaInPixel')
hmaRgiPixelDir = path.join(stepsDir, 'Step6_HMA_rgi60d_Pixel')
hmaPixelStaDir = path.join(stepsDir, 'Step7_HMA_Pixel_Sta')

hmaRigPixelPath = path.join(hmaRgiPixelDir, 'HMA_rgi60d_ModisPixelD_' + cellsize + '.shp')

hmaPixelStaName = 'HMA_ModisPixel_{0}_sta_{1}.shp'.format(cellsize, topoType)
hmaPixelStaPath = path.join(hmaPixelStaDir, hmaPixelStaName)
if not arcpy.Exists(hmaPixelStaPath):
    arcpy.CopyFeatures_management(hmaRigPixelPath, hmaPixelStaPath)

for staType in staTypeList:
    topoStaTableName = hmaPixelStaName.replace('.shp', '_' + staType + '.dbf')
    topoStaTablePath = path.join(hmaPixelStaDir, topoStaTableName)
    if not arcpy.Exists(topoStaTablePath):
        arcpy.sa.ZonalStatisticsAsTable(hmaPixelStaPath, 'FID', srtmPath, topoStaTablePath, 'DATA',
                                        staType.upper())
    staType = staType.rstrip('imum').upper()
    fieldList = [field.name for field in arcpy.ListFields(hmaPixelStaPath)]
    if staType not in fieldList:
        arcpy.JoinField_management(hmaPixelStaPath, 'FID', topoStaTablePath, 'FID_', [staType])
