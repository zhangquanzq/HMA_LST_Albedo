# coding=utf-8
# 检查HMA地区的SRTM 1arc数据是否正常存在。
import arcpy
import glob
import os
import zipfile

fishnetHmaPath = 'F:\\HMA_LST_Albedo\\Data\\Basic\\Feature\\SRTM_1arc\\Fishnet_SRTM_1arc_HMA.shp'
srtmCnTileZipFolderPath = 'F:\\HMA_LST_Albedo\\Data\\Basic\\Raster\\SRTM1arc_Zips\\'

restSrtmHmaTileZipPath = 'F:\\HMA_LST_Albedo\\ss.txt'

srtmHmaTileIdList = [record.getValue('TileID') for record in arcpy.SearchCursor(fishnetHmaPath, fields='TileID')]

srtmCnTileZipPathList = glob.glob(os.path.join(srtmCnTileZipFolderPath, '*.zip'))
srtmCnTileZipList = [os.path.basename(srtmCnZipFile) for srtmCnZipFile in srtmCnTileZipPathList]

validSrtmHmaTileZipList, restSrtmHmaTileZipList = [], []
for tileId in srtmHmaTileIdList:
    srtmTileZipFileName = tileId + '.SRTMGL1.hgt.zip'
    if srtmTileZipFileName in srtmCnTileZipList:
        validSrtmHmaTileZipList.append(srtmTileZipFileName)
    else:
        restSrtmHmaTileZipList.append(srtmTileZipFileName + '\n')

valideSrtmHmaTileZipPathList = []
for srtmHmaZipFile in validSrtmHmaTileZipList:
    valideSrtmHmaTileZipPathList.append(os.path.join(srtmCnTileZipFolderPath, srtmHmaZipFile))

with open(restSrtmHmaTileZipPath, 'w') as f:
    f.writelines(restSrtmHmaTileZipList)

for srtmTileZipPath in valideSrtmHmaTileZipPathList:
    srtmHgtFile = zipfile.ZipFile(srtmTileZipPath)
    srtmHgtFile.extractall('F:\\HMA_LST_Albedo\\Data\\Basic\\Raster\\SRTM1arc_Tiles_HMA\\')
