# coding=utf-8
import arcpy
import os
from os import path

# 确定分辨率的标记. 1表示500m, 2表示1km.
flg1 = 2

cellsize = ['500m', '1km'][flg1 - 1]

# 路径.
rootDir = r'G:\HMA_LST_Albedo\Data'
refPath = [path.join(rootDir, r'MOD10A1\MOD10A1_2_MosaicHMA_TIF\MOD10A1_2001XXX_TIF',
                     'MOD10A1.A2001001.Albedo.tif'),
           path.join(rootDir, r'MOD11A1\MOD11A1_2_MosaicCN_TIF\MOD11A1_2001XXX_TIF',
                     'MOD11A1.A2001001.LST_Day.tif')][flg1 - 1]
stepsDir = path.join(rootDir, 'GlacierAreaInPixel')
pctFeatureDir = path.join(stepsDir, 'Step8_HMA_Pixel_Sta_Percent')
pctRasterDir = path.join(stepsDir, 'Step9_HMA_Pixel_Percent_Raster')
if not path.exists(pctRasterDir):
    os.mkdir(pctRasterDir)

# 矢量转栅格.
refCellsize = arcpy.Raster(refPath).meanCellHeight
arcpy.env.workspace = pctFeatureDir
arcpy.env.snapRaster = refPath
arcpy.env.extent = arcpy.Raster(refPath).extent

featureList = arcpy.ListFeatureClasses('HMA*{0}*elev*.shp'.format(cellsize))
for feature in featureList:
    pixelRasterName = feature.replace('sta_elev', 'rgi60').replace('.shp', '.tif')
    pixelRasterPath = path.join(pctRasterDir, pixelRasterName)
    if not arcpy.Exists(pixelRasterPath):
        arcpy.PolygonToRaster_conversion(feature, 'Id', pixelRasterPath, cellsize=refCellsize)
        print(pixelRasterName)

    pixelFidRasterName = feature.replace('sta_elev', 'rgi60_FID').replace('.shp', '.tif')
    pixelFidRasterPath = path.join(pctRasterDir, pixelFidRasterName)
    if not arcpy.Exists(pixelFidRasterPath):
        arcpy.PolygonToRaster_conversion(feature, 'FID', pixelFidRasterPath, cellsize=refCellsize)
        print(pixelFidRasterName)
