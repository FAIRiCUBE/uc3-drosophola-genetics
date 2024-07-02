import csv
#from code.module_crs_converter import trans4mEPSG
from code.functions import trans4mEPSG


class Coverage(object):
    def __init__(self,data):
        self._data=data
        self._ID= data[0]
        self._CRS=data[1]
        self._minlat=data[2]
        self._maxlar= data[3]
        self._minlong=data[4]
        self._maxlong=data[5]
        self._fromtime=data[6] 
        self._totime=data[7]
        lnames=[]
        for lna in range(1,len(data)):
            lnames.append(data[lna][0])
        self.layers=lnames
    def getBoundary(self):
        minlats=[]
        maxlats=[]
        minlongs=[]
        maxlongs=[]  
        for minl in range(1,len(self._data)):
            crs=self._data[minl][1]
            if crs=="EPSG/0/4326":
                min_y, min_x = trans4mEPSG("EPSG:4326","EPSG:3035", float(self._data[minl][4]), float(self._data[minl][2]))
                max_y, max_x = trans4mEPSG("EPSG:4326","EPSG:3035",float(self._data[minl][5]), float(self._data[minl][3]))
                minlats.append(min_x)
                maxlats.append(max_x)
                minlongs.append(min_y)
                maxlongs.append(max_y)
            else:
                minlats.append(float(self._data[minl][2]))
                maxlats.append(float(self._data[minl][3]))
                minlongs.append(float(self._data[minl][4]))
                maxlongs.append(float(self._data[minl][5]))
        self.minlat=min([x for x in minlats if x != float('inf')])
        self.maxlat= max([x for x in maxlats if x != float('inf')])
        self.minlong=min([x for x in minlongs if x != float('inf')])
        self.maxlong=max([x for x in maxlongs if x != float('inf')])
    def getSamples(self,path):
        filtered_data={}
        with open(path, 'r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                lat=row["lat"]
                long=row["long"]
                if lat == 'NA' and long == 'NA':
                    continue
                #else:
                    #long,lat=trans4mEPSG("EPSG:4326","EPSG:3035",float(long),float(lat))
                if float(lat) > self.minlat and float(lat) < self.maxlat and float(long) > self.minlong and float(long) < self.maxlong:
                    sampleinfo=(row["lat"],row["long"])
                    filtered_data[row["sampleId"]] = sampleinfo
                    #filtered_data.append(sampleinfo)
        self.samples=filtered_data
        