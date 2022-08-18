# coding=utf-8
import arcpy
import numpy as np
from os import path

# 标记和预设参数.
# 指定像元分辨率(数据类型)的标记. 1表示1km(MODIS LST), 2表示500m(MODIS Albedo).
flg1 = 1

cellsize = ['1km', '500m'][flg1 - 1]
areaTotal = ['858635', '214658.673273'][flg1 - 1]

regionList = ['Aertaishan', 'Kunlunshan', 'Qilianshan', 'Tanggulashan', 'Tianshan', 'Ximalaya',
              'Zangnan']
topoList = ['elev', 'slp']
arcpy.env.qualifiedFieldNames = False

# 路径.
dataDir = path.join(r'F:\HMA_LST_Albedo\Data')
glacierInPixelDir = path.join(dataDir, 'GlacierAreaInPixel')

hmaExtentDir = path.join(glacierInPixelDir, 'Step1_HMA_Extent_Subregions')
hmaRgiRegionDir = path.join(glacierInPixelDir, 'Step2_HMA_rgi60_Subregions')
hmaRgiDRegionDir = path.join(glacierInPixelDir, 'Step3_HMA_rgi60d')
hmaFishnetDir = path.join(glacierInPixelDir, 'Step4_HMA_Fishnet_Subregions')
hmaPixelRegionDir = path.join(glacierInPixelDir, 'Step5_HMA_Pixel_Subregions')
hmaRgiPixelRegionDir = path.join(glacierInPixelDir, 'Step6_HMA_Pixel_rgi60d_Subregions')
hmaPixelDir = path.join(glacierInPixelDir, 'Step7_HMA_Pixel')
hmaRgiPixelDir = path.join(glacierInPixelDir, 'Step8_HMA_rgi60d_Pixel')
hmaPixelStaDir = path.join(glacierInPixelDir, 'Step9_HMA_Pixel_Sta')
hmaPixelStaPctDir = path.join(glacierInPixelDir, 'Step10_HMA_Pixel_Sta_Percent')
hmaPixelPctRasterDir = path.join(glacierInPixelDir, 'Step11_HMA_Pixel_Percent_Raster')

featurePath = path.join(dataDir, r'Basic\Feature\00_rgi60_O2Regions.shp')

refModisDir = [r'MOD11A1\MOD11A1_2_MosaicCN_TIF\MOD11A1_2001XXX_TIF',
               r'MOD10A1\MOD10A1_2_MosaicHMA_TIF\MOD10A1_2001XXX_TIF'][flg1 - 1]
refModisName = ['MOD11A1.A2001001.LST_Day.tif', 'MOD10A1.A2001001.Albedo.tif'][flg1 - 1]
refModisPath = path.join(dataDir, refModisDir, refModisName)

refModisRas = arcpy.Raster(refModisPath)
refExtent = refModisRas.extent
cellHeight, cellWidth = refModisRas.meanCellHeight, refModisRas.meanCellWidth

# 数据处理过程.
hmaSubregionPath = path.join(hmaExtentDir, 'HMA_Subregions.shp')
hmaRgiPath = path.join(hmaRgiRegionDir, 'HMA_rgi60.shp')
for region in regionList:
    # Step1, 提取HMA分区矢量中每个分区多边形为单独的矢量文件.
    regionExtentPath = path.join(hmaExtentDir, 'HMA_Extent_{0}.shp'.format(region))
    if not arcpy.Exists(regionExtentPath):
        print('导出{0}区的边界矢量.'.format(region))
        arcpy.Select_analysis(hmaSubregionPath, regionExtentPath, '"Name" = \'{0}\''.format(region))

    # Step2, 按HMA每个分区的范围提取各分区的冰川边界矢量文件.
    hmaRgiRegionPath = path.join(hmaRgiRegionDir, 'HMA_rgi60_{0}.shp'.format(region))
    if not arcpy.Exists(hmaRgiRegionPath):
        print('提取{0}区的Rgi60边界.'.format(region))
        arcpy.Clip_analysis(hmaRgiPath, regionExtentPath, hmaRgiRegionPath)

    # Step3a, 消除邻接冰川的边界线, 防止提取冰川覆盖像元时, 漏选边界线上的像元.
    hmaRgiDRegionPath = path.join(hmaRgiDRegionDir, 'HMA_rgi60d_{0}.shp'.format(region))
    if not arcpy.Exists(hmaRgiDRegionPath):
        print('消除{0}区Rgi60的内边界.'.format(region))
        arcpy.Dissolve_management(hmaRgiRegionPath, hmaRgiDRegionPath, 'EndDate', '', 'SINGLE_PART')

# Step3b, 合并消除内边界的冰川边界范围矢量.
hmaRgiDPath = path.join(hmaRgiDRegionDir, 'HMA_rgi60d.shp')
if not arcpy.Exists(hmaRgiDPath):
    print('拼接HMA所有分区的消除过内边界的Rgi60边界矢量.')
    arcpy.env.workspace = hmaRgiDRegionDir
    hmaRgiDRegionList = arcpy.ListFeatureClasses('HMA_rgi60d_*.shp')
    arcpy.Merge_management(hmaRgiDRegionList, hmaRgiDPath)

for region in regionList:
    regionCellsize = '{0}_{1}'.format(region, cellsize)
    # Step4, 为每个分区多边形创建对应的MODIS像元Fishnet矢量文件.
    regionExtentPath = path.join(hmaExtentDir, 'HMA_Extent_{0}.shp'.format(region))
    regionFishnetPath = path.join(hmaFishnetDir, 'HMA_Fishnet_{0}.shp'.format(regionCellsize))
    if not arcpy.Exists(regionFishnetPath):
        print('创建{0}的Fishnet.'.format(regionCellsize))
        featureExtent = arcpy.da.SearchCursor(regionExtentPath, ['SHAPE@']).next()[0].extent
        xMinFishnet = np.arange(refExtent.XMin, featureExtent.XMin, cellWidth)[-1]
        xMaxFishnet = np.arange(refExtent.XMax, featureExtent.XMax, -cellWidth)[-1]
        yMinFishnet = np.arange(refExtent.YMin, featureExtent.YMin, cellHeight)[-1]
        yMaxFishnet = np.arange(refExtent.YMax, featureExtent.YMax, -cellHeight)[-1]
        lowerLeft = '{0} {1}'.format(xMinFishnet, yMinFishnet)
        upperLeft = '{0} {1}'.format(xMinFishnet, yMaxFishnet)
        upperRight = '{0} {1}'.format(xMaxFishnet, yMaxFishnet)
        arcpy.env.outputCoordinateSystem = arcpy.Describe(refModisPath).SpatialReference
        arcpy.CreateFishnet_management(regionFishnetPath, lowerLeft, upperLeft, cellWidth,
                                       cellHeight, '#', '#', upperRight, 'NO_LABELS', '#',
                                       'POLYGON')
        arcpy.ClearEnvironment('outputCoordinateSystem')

    # Step5, 提取与各HMA分区矢量范围重叠或相交的MODIS像元Fishnet格网.
    hmaRgiDRegionPath = path.join(hmaRgiDRegionDir, 'HMA_rgi60d_{0}.shp'.format(region))
    hmaPixelRegionName = 'HMA_ModisPixel_{0}.shp'.format(regionCellsize)
    hmaPixelRegionPath = path.join(hmaPixelRegionDir, hmaPixelRegionName)
    if not arcpy.Exists(hmaPixelRegionPath):
        print('提取{0}的MODIS像元矢量.'.format(regionCellsize))
        layerName = 'pixelLayer_{0}'.format(regionCellsize)
        arcpy.MakeFeatureLayer_management(regionFishnetPath, layerName)
        arcpy.SelectLayerByLocation_management(layerName, 'INTERSECT', hmaRgiDRegionPath)
        arcpy.CopyFeatures_management(layerName, hmaPixelRegionPath)

    # Step6a, 叠置Pixel和Rgi60矢量图层, 并计算叠置后各多边形的面积.
    hmaRgiRegionPath = path.join(hmaRgiRegionDir, 'HMA_rgi60_{0}.shp'.format(region))
    pixelRgiPath = path.join(hmaRgiPixelRegionDir, 'HMA_Pixel_rgi60_{0}.shp'.format(regionCellsize))
    if not arcpy.Exists(pixelRgiPath):
        print('叠置{0}的像元图层和Rgi60图层.'.format(regionCellsize))
        arcpy.Identity_analysis(hmaPixelRegionPath, hmaRgiRegionPath, pixelRgiPath)
    if 'Area_1' not in [field.name for field in arcpy.ListFields(pixelRgiPath)]:
        arcpy.AddField_management(pixelRgiPath, 'Area_1', 'DOUBLE')
        arcpy.CalculateField_management(pixelRgiPath, 'Area_1', '!shape.area!', 'PYTHON_9.3')

    # Step6b, 获取每个有冰川Pixel的编号, 及其范围内面积最大多边形所属冰川的编号及面积, 并保存为CSV文件.
    csvPath = pixelRgiPath.replace('.shp', '_selectRecord.csv')
    if not path.exists(csvPath):
        print('筛选每个{0}像元内面积最大的冰川信息.'.format(regionCellsize))
        fieldList1 = ['FID_HMA_Mo', 'FID_HMA_rg', 'RGIId', 'Area_1']
        pixelRgiTable = arcpy.da.TableToNumPyArray(pixelRgiPath, fieldList1)
        fidHmaMoList = pixelRgiTable['FID_HMA_Mo']
        fidHmaRgList = pixelRgiTable['FID_HMA_rg']
        rgiidList = pixelRgiTable['RGIId']
        area1List = pixelRgiTable['Area_1']

        selectedRecordList = ['mFID,RGIId,maxArea\n']
        for fidHmaMo in np.unique(fidHmaMoList):
            fidHmaMoIndex = np.where(fidHmaMo == fidHmaMoList)
            fidHmaRgVector = fidHmaRgList[fidHmaMoIndex]

            # 当Pixel范围内没有冰川时(冰川矢量边界与像元矢量边界邻接但不相交), 跳过.
            if np.size(fidHmaMoIndex) == 1 and -1 in fidHmaRgVector:
                continue
            rgiidVector = rgiidList[fidHmaMoIndex]
            area1Vector = area1List[fidHmaMoIndex]

            areaMax = max(np.delete(area1Vector, np.where(fidHmaRgVector == -1)))
            rgiid = rgiidVector[np.where(area1Vector == areaMax)[0][0]]  # 若最大面积不止1个, 取第1个.
            selectedRecordList.append('{0},{1},{2}\n'.format(fidHmaMo, rgiid, areaMax))
        with open(csvPath, 'w') as csvf:
            csvf.writelines(selectedRecordList)

    # Step7a, 将CSV文件中的属性与HMA冰川像元矢量属性表链接.
    maxAreaShpPath = path.join(hmaPixelDir, 'HMA_ModisPixel_{0}_maxArea.shp'.format(regionCellsize))
    dbfPath = csvPath.replace('.csv', '.dbf')
    if not arcpy.Exists(maxAreaShpPath):
        layerName = regionCellsize
        print('链接面积最大冰川信息与对应的{0} MODIS像元.'.format(layerName))
        arcpy.MakeFeatureLayer_management(hmaPixelRegionPath, layerName)
        arcpy.CopyRows_management(csvPath, dbfPath)
        arcpy.AddJoin_management(layerName, 'FID', dbfPath, 'mFID')
        arcpy.CopyFeatures_management(layerName, maxAreaShpPath)
        arcpy.DeleteField_management(maxAreaShpPath, 'OID_')
        arcpy.Delete_management(dbfPath)

# Step7b: 拼接HMA分区文件, 空间链接Rgi分区属性字段.
pixelMergePath = path.join(hmaPixelDir, 'HMA_ModisPixel_{0}_Merge.shp'.format(cellsize))
if not arcpy.Exists(pixelMergePath):
    print('拼接HMA所有分区的{0}像元图层.'.format(cellsize))
    arcpy.env.workspace = hmaPixelDir
    maxAreaShpList = arcpy.ListFeatureClasses('*{0}_maxArea.shp'.format(cellsize))
    arcpy.Merge_management(maxAreaShpList, pixelMergePath)

pixelRgiRegionName = 'HMA_ModisPixel_{0}_RgiRegion.shp'.format(cellsize)
sJoinRgiRegionPath = path.join(hmaPixelDir, pixelRgiRegionName)
if not arcpy.Exists(sJoinRgiRegionPath):
    print('链接Rgi编号到HMA{0}像元图层.'.format(cellsize))
    arcpy.SpatialJoin_analysis(pixelMergePath, featurePath, sJoinRgiRegionPath)

# Step8, 将消除了内边界的HMA冰川范围矢量与Pixel矢量叠置, 并合并每个Pixel范围内的多边形.
hmaRgiDPixelPath = path.join(hmaRgiPixelDir, 'HMA_rgi60d_ModisPixel_{0}.shp'.format(cellsize))
if not arcpy.Exists(hmaRgiDPixelPath):
    arcpy.Identity_analysis(hmaRgiDPath, sJoinRgiRegionPath, hmaRgiDPixelPath)

hmaRgiDPixelDPath = path.join(hmaRgiPixelDir, 'HMA_rgi60d_ModisPixelD_{0}.shp'.format(cellsize))
if not arcpy.Exists(hmaRgiDPixelDPath):
    arcpy.Dissolve_management(hmaRgiDPixelPath, hmaRgiDPixelDPath, 'FID_HMA_Mo', '', 'MULTI_PART')
fieldList2 = [field.name for field in arcpy.ListFields(hmaRgiDPixelDPath)]
if 'Area' not in fieldList2:
    arcpy.AddField_management(hmaRgiDPixelDPath, 'Area', 'DOUBLE')
    arcpy.CalculateField_management(hmaRgiDPixelDPath, 'Area', '!shape.area!', 'PYTHON_9.3')
if 'Area_Pct' not in fieldList2:
    arcpy.AddField_management(hmaRgiDPixelDPath, 'Area_Pct', 'DOUBLE')
    arcpy.CalculateField_management(hmaRgiDPixelDPath, 'Area_Pct', '!Area!/{0}'.format(areaTotal),
                                    'PYTHON_9.3')

# Step9: 统计冰川像元高程, 坡度值, 链接有用字段, 删除无用字段.
validFieldList = ['FID', 'Shape', 'Id', 'mFID', 'RGIId', 'maxArea', 'FULL_NAME', 'RGI_CODE',
                  'WGMS_CODE', 'FID_HMA_Mo', 'Area', 'Area_Pct', 'FID_HMA__1', 'COUNT', 'MIN',
                  'MAX', 'RANGE', 'MEAN', 'STD']
rgiDPixelDPath = path.join(hmaRgiPixelDir, 'HMA_rgi60d_ModisPixelD_{0}.shp'.format(cellsize))
for topo in topoList:
    # 统计冰川像元高程, 坡度值.
    topoPath = path.join(dataDir, r'Basic\Raster\HMA_SRTM_3arc_{0}.tif'.format(topo))
    pixelStaTableName = 'HMA_ModisPixel_{0}_Sta_{1}_All.dbf'.format(cellsize, topo)
    pixelStaTablePath = path.join(hmaPixelStaDir, pixelStaTableName)
    if not arcpy.Exists(pixelStaTablePath):
        print('分区统计{0}像元的{1}.'.format(cellsize, topo))
        arcpy.sa.ZonalStatisticsAsTable(rgiDPixelDPath, 'FID_HMA_Mo', topoPath, pixelStaTablePath)

    # 链接有用字段(每次只能连接一个表格).
    pixelStaShpPath = pixelStaTablePath.replace('_All.dbf', '.shp')
    if not arcpy.Exists(pixelStaShpPath):
        tempShpPath = pixelStaShpPath.replace('.shp', '_tmp.shp')
        if arcpy.Exists(tempShpPath):
            arcpy.Delete_management(tempShpPath)

        layerName1 = 'layer7{0}'.format(topo)
        arcpy.MakeFeatureLayer_management(sJoinRgiRegionPath, layerName1)
        arcpy.AddJoin_management(layerName1, 'FID', rgiDPixelDPath, 'FID_HMA_Mo')
        arcpy.CopyFeatures_management(layerName1, tempShpPath)

        layerName2 = 'layer7b{0}'.format(topo)
        arcpy.MakeFeatureLayer_management(tempShpPath, layerName2)
        arcpy.AddJoin_management(layerName2, 'FID', pixelStaTablePath, 'FID_HMA_Mo')
        arcpy.CopyFeatures_management(layerName2, pixelStaShpPath)

        arcpy.Delete_management(tempShpPath)

    # 删除无用字段.
    staShpFieldList = [field.name for field in arcpy.ListFields(pixelStaShpPath)]
    uselessFieldList = list(set(staShpFieldList) - set(validFieldList))
    if len(uselessFieldList) > 0:
        arcpy.DeleteField_management(pixelStaShpPath, uselessFieldList)

    # 导出像元内面积比例大于80%的像元.
    hmaPixelStaPctName = 'HMA_ModisPixel_{0}_sta_{1}_80percent.shp'.format(cellsize, topo)
    hmaPixelStaPctPath = path.join(hmaPixelStaPctDir, hmaPixelStaPctName)
    if not arcpy.Exists(hmaPixelStaPctPath):
        print('导出面积百分比80%以上的{0} {1}统计像元.'.format(cellsize, topo))
        arcpy.Select_analysis(pixelStaShpPath, hmaPixelStaPctPath, '"Area_Pct" >= 0.8')

    # 矢量转栅格.
    hmaPixelPctRasterName = 'HMA_ModisPixel_{0}_rgi60_FID_80percent.tif'.format(cellsize)
    hmaPixelPctRasterPath = path.join(hmaPixelPctRasterDir, hmaPixelPctRasterName)
    if not arcpy.Exists(hmaPixelPctRasterPath):
        arcpy.env.extent = refModisPath
        arcpy.env.SnapRaster = refModisPath
        arcpy.FeatureToRaster_conversion(hmaPixelStaPctPath, 'FID', hmaPixelPctRasterPath,
                                         cell_size=cellHeight)
