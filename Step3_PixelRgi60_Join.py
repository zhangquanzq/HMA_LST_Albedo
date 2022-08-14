# coding=utf-8
import arcpy
import numpy as np
from os import path

# 标记和预设参数.
# 指定像元分辨率(数据类型)的标记. 1表示1km(MODIS LST), 2表示500m(MODIS Albedo).
flg1 = 1

arcpy.env.qualifiedFieldNames = False

cellsize = ['1km', '500m'][flg1 - 1]

regionList = ['Aertaishan', 'Kunlunshan', 'Qilianshan', 'Tanggulashan', 'Tianshan', 'Ximalaya',
              'Zangnan']
topoList = ['elev', 'slp']

# 路径.
dataDir = path.join(r'F:\HMA_LST_Albedo\Data')
glacierInPixelDir = path.join(dataDir, 'GlacierAreaInPixel')

hmaRgiRegionDir = path.join(glacierInPixelDir, 'Step3_HMA_rgi60_Subregions')
hmaPixelRegionDir = path.join(glacierInPixelDir, 'Step5_HMA_Pixel_Subregions')
hmaRgiPixelRegionDir = path.join(glacierInPixelDir, 'Step6_HMA_rgi60d_Pixel_Subregions')
hmaRgiPixelDir = path.join(glacierInPixelDir, 'Step6_HMA_rgi60d_Pixel')
hmaPixelStaDir = path.join(glacierInPixelDir, 'Step7_HMA_Pixel_Sta')

featurePath = path.join(dataDir, r'Basic\Feature\00_rgi60_O2Regions.shp')

# 处理步骤.
# Step1, 提取HMA分区矢量中每个分区多边形为单独的矢量文件.

# Step2, 为每个分区多边形创建对应的MODIS像元Fishnet矢量文件.

# Step3, 按HMA每个分区的范围提取各分区的冰川边界矢量文件.

# Step4, 消除邻接冰川的边界线, 防止提取冰川覆盖像元时, 漏选边界线上的像元.

# Step5, 提取与各HMA分区矢量范围重叠或相交的MODIS像元Fishnet格网.

# Step6, 将


for region in regionList:
    regionCellsize = '{0}_{1}'.format(region, cellsize)
    rgiPath = path.join(hmaRgiRegionDir, 'HMA_rgi60_{0}.shp'.format(region))
    pixelPath = path.join(hmaPixelRegionDir, 'HMA_ModisPixel_{0}.shp'.format(regionCellsize))

    pixelRgiPath = path.join(hmaRgiPixelRegionDir, 'HMA_Pixel_rgi60_{0}.shp'.format(regionCellsize))
    pixelJoinPath = pixelPath.replace('.shp', '_maxArea.shp')
    csvPath = pixelRgiPath.replace('.shp', '_selectRecord.csv')
    dbfPath = csvPath.replace('.csv', '.dbf')

    # Step6
    # 叠置分析Pixel和Rgi60图层数据, 并计算叠置后各多边形的面积.
    if not arcpy.Exists(pixelRgiPath):
        print(u'叠置Pixel和Rgi60图层: {0}'.format(regionCellsize))
        arcpy.Identity_analysis(pixelPath, rgiPath, pixelRgiPath)
    if 'Area_1' not in [field.name for field in arcpy.ListFields(pixelRgiPath)]:
        arcpy.AddField_management(pixelRgiPath, 'Area_1', 'DOUBLE')
        arcpy.CalculateField_management(pixelRgiPath, 'Area_1', '!shape.area!', 'PYTHON_9.3')

    # 获取每个Pixel的编号, 及其范围内面积最大多边形所属冰川的编号及面积, 并保存为CSV文件.
    if not path.exists(csvPath):
        print(u'获取像元内面积最大冰川编号: {0}'.format(regionCellsize))
        # 置后矢量数据属性表有用字段记录.
        sFieldList = ['FID_HMA_Mo', 'FID_HMA_rg', 'RGIId', 'Area_1']
        pixelRgiArray = arcpy.da.TableToNumPyArray(pixelRgiPath, sFieldList)
        fidHmaMoList = pixelRgiArray['FID_HMA_Mo']
        fidHmaRgList = pixelRgiArray['FID_HMA_rg']
        rgiidList = pixelRgiArray['RGIId']
        area1List = pixelRgiArray['Area_1']

        selectedRecordList = ['mFID,RGIId,maxArea\n']
        for fidHmaMo in np.unique(fidHmaMoList):
            fidHmaMoIndex = np.where(fidHmaMo == fidHmaMoList)
            fidHmaRgVector = fidHmaRgList[fidHmaMoIndex]

            # 当Pixel范围内没有冰川时(冰川矢量边界与像元矢量边界相接但不相交), 跳过.
            if np.size(fidHmaMoIndex) == 1 and -1 in fidHmaRgVector:
                continue
            rgiidVector = rgiidList[fidHmaMoIndex]
            area1Vector = area1List[fidHmaMoIndex]

            areaMax = max(np.delete(area1Vector, np.where(fidHmaRgVector == -1)))
            rgiid = rgiidVector[np.where(area1Vector == areaMax)[0][0]]  # 若最大面积不止1个, 取第1个.
            selectedRecordList.append('{0},{1},{2}\n'.format(fidHmaMo, rgiid, areaMax))
        with open(csvPath, 'w') as csvf:
            csvf.writelines(selectedRecordList)

    if not arcpy.Exists(pixelJoinPath):
        layerName = regionCellsize
        print(u'链接面积最大冰川编号到Pixel图层属性表: {0}'.format(layerName))
        arcpy.MakeFeatureLayer_management(pixelPath, layerName)
        arcpy.CopyRows_management(csvPath, dbfPath)
        arcpy.AddJoin_management(layerName, 'FID', dbfPath, 'mFID')
        arcpy.CopyFeatures_management(layerName, pixelJoinPath)
        arcpy.DeleteField_management(pixelJoinPath, 'OID_')
        arcpy.Delete_management(dbfPath)

# Step5: 拼接HMA分区文件, 空间链接Rgi分区属性字段.
mergeMaxAreaPath = path.join(hmaPixelRegionDir, 'HMA_ModisPixel_{0}_Merge.shp'.format(cellsize))
if not arcpy.Exists(mergeMaxAreaPath):
    print(u'拼接HMA各分区Pixel图层: {0}'.format(cellsize))
    arcpy.env.workspace = hmaPixelRegionDir
    maxAreaShpList = arcpy.ListFeatureClasses('*{0}_maxArea.shp'.format(cellsize))
    arcpy.Merge_management(maxAreaShpList, mergeMaxAreaPath)

sJoinRgiRegionName = 'HMA_ModisPixel_{0}_RgiRegion.shp'.format(cellsize)
sJoinRgiRegionPath = path.join(hmaPixelRegionDir, sJoinRgiRegionName)
if not arcpy.Exists(sJoinRgiRegionPath):
    arcpy.SpatialJoin_analysis(mergeMaxAreaPath, featurePath, sJoinRgiRegionPath)

# Step7: 统计冰川像元高程, 坡度值, 链接有用字段, 删除无用字段.
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
        print('统计{0}像元的{1}'.format(cellsize, topo))
        arcpy.sa.ZonalStatisticsAsTable(rgiDPixelDPath, 'FID_HMA_Mo', topoPath, pixelStaTablePath)

    # 链接有用字段.
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
