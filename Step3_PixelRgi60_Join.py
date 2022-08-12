# coding=utf-8
import arcpy
import numpy as np
from os import path

# 标记和预设参数.
# 指定像元分辨率(数据类型)的标记. 1表示1km(MODIS LST), 2表示500m(MODIS Albedo).
flg1 = 2

arcpy.env.qualifiedFieldNames = False

cellsize = ['1km', '500m'][flg1 - 1]

regionList = ['Aertaishan', 'Kunlunshan', 'Qilianshan', 'Tanggulashan', 'Tianshan', 'Ximalaya',
              'Zangnan']

# 路径.
rootDir = r'F:\HMA_LST_Albedo'
glacierInPixelDir = path.join(rootDir, r'Data\GlacierAreaInPixel')

hmaRgiRegionDir = path.join(glacierInPixelDir, 'Step3_HMA_rgi60_Subregions')
hmaPixelRegionDir = path.join(glacierInPixelDir, 'Step5_HMA_Pixel_Subregions')
hmaRgiPixelRegionDir = path.join(glacierInPixelDir, 'Step6_HMA_rgi60d_Pixel_Subregions')

# 处理每个HMA分区.
for region in regionList:
    regionCellsize = '{0}_{1}'.format(region, cellsize)
    rgiPath = path.join(hmaRgiRegionDir, 'HMA_rgi60_{0}.shp'.format(region))
    pixelPath = path.join(hmaPixelRegionDir, 'HMA_ModisPixel_{0}.shp'.format(regionCellsize))

    pixelRgiPath = path.join(hmaRgiPixelRegionDir, 'HMA_Pixel_rgi60_{0}.shp'.format(regionCellsize))
    pixelJoinPath = pixelPath.replace('.shp', '_Join.shp')
    csvPath = pixelRgiPath.replace('.shp', '_SelectRecord.csv')
    dbfPath = csvPath.replace('.csv', '.dbf')

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
        print(u'链接面积最大冰川编号到Pixel图层属性表: {0}'.format(regionCellsize))
        layerName = regionCellsize
        arcpy.MakeFeatureLayer_management(pixelPath, regionCellsize)
        arcpy.CopyRows_management(csvPath, dbfPath)
        arcpy.AddJoin_management(regionCellsize, 'FID', dbfPath, 'mFID')
        arcpy.CopyFeatures_management(regionCellsize, pixelJoinPath)
        arcpy.DeleteField_management(pixelJoinPath, 'OID_')
        arcpy.Delete_management(dbfPath)

# 拼接各HMA分区的Join文件.
joinMergePath = path.join(hmaPixelRegionDir, 'HMA_ModisPixel_{0}_Join.shp'.format(cellsize))
if not arcpy.Exists(joinMergePath):
    print(u'拼接HMA各分区Pixel图层: {0}'.format(cellsize))
    arcpy.env.workspace = hmaPixelRegionDir
    joinNameList = arcpy.ListFeatureClasses('*Join.shp')
    arcpy.Merge_management(joinNameList, joinMergePath)
