# coding=utf-8
import arcpy
import numpy as np
from os import path

# 路径.
rootDir = r'E:\HMA_LST_Albedo'
dataDir = path.join(rootDir, 'Data')
stepsDir = path.join(dataDir, 'GlacierAreaInPixel')

testDir = path.join(stepsDir, 'test')
rgi60Path = path.join(testDir, 'HMA_rgi60_Aertaishan2.shp')
pixelPath = path.join(testDir, 'HMA_ModisPixel_Aertaishan_1km.shp')

pixelRgi60Path = path.join(testDir, 'HMA_Pixel_rgi60_Aertaishan.shp')
csvFilePath = path.join(testDir, 'Export_Output2.dbf')

# 叠置分析Pixel和Rgi60图层数据, 并计算叠置后各多边形的面积.
if not arcpy.Exists(pixelRgi60Path):
    arcpy.Identity_analysis(pixelPath, rgi60Path, pixelRgi60Path)
allFieldList = [field.name for field in arcpy.ListFields(pixelRgi60Path)]
if 'Area_1' not in allFieldList:
    arcpy.AddField_management(pixelRgi60Path, 'Area_1', 'DOUBLE')
    arcpy.CalculateField_management(pixelRgi60Path, 'Area_1', '!shape.area!', 'PYTHON_9.3')

# 提取叠置后矢量数据属性表有用字段.
selectedFieldList = ['FID_HMA_Mo', 'FID_HMA_rg', 'RGIId', 'Area_1']
pixelRgiArray = arcpy.da.TableToNumPyArray(pixelRgi60Path, selectedFieldList)
fidHmaMoList = pixelRgiArray['FID_HMA_Mo']
fidHmaRgList = pixelRgiArray['FID_HMA_rg']
rgiidList = pixelRgiArray['RGIId']
area1List = pixelRgiArray['Area_1']

# 获取每个Pixel的编号,及其范围内面积最大多边形所属冰川的编号及面积, 并保存为CSV文件.
if not path.exists(csvFilePath):
    selectedRecordList = ['mFID,RGIId,maxArea\n']
    for fidHmaMo in np.unique(fidHmaMoList):
        fidHmaMoIndex = np.where(fidHmaMo == fidHmaMoList)

        fidHmaRgVector = fidHmaRgList[fidHmaMoIndex]
        rgiidVector = rgiidList[fidHmaMoIndex]
        area1Vector = area1List[fidHmaMoIndex]

        areaMax = max(np.delete(area1Vector, np.where(fidHmaRgVector == -1)))
        rgiid = rgiidVector[np.where(area1Vector == areaMax)[0][0]]  # 若最大面积多边形不止1个, 取第1个.
        selectedRecordList.append('{0},{1},{2}\n'.format(fidHmaMo, rgiid, areaMax))
    with open(csvFilePath, 'w') as csvf:
        csvf.writelines(selectedRecordList)

arcpy.AddJoin_management(pixelPath, 'FID', csvFilePath, 'OID')
arcpy.Copy_management(pixelPath, 'd:\\ddd.shp')

d = 1
