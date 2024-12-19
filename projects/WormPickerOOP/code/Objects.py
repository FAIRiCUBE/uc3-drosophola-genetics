import csv
#from code.module_crs_converter import trans4mEPSG
from code.functions import trans4mEPSG
import logging


class Coverage(object):
    def __init__(self,data):
        self._data=data
        self._ID= data[0]
        self._CRS=data[1]
        self._minlat=data[2]
        self._maxlat= data[3]
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
                #min_y, min_x = trans4mEPSG("EPSG:4326","EPSG:3035", float(self._data[minl][4]), float(self._data[minl][2]))
                #max_y, max_x = trans4mEPSG("EPSG:4326","EPSG:3035",float(self._data[minl][5]), float(self._data[minl][3]))
                minlats.append(float(self._data[minl][2]))
                maxlats.append(float(self._data[minl][3]))
                minlongs.append(float(self._data[minl][4]))
                maxlongs.append(float(self._data[minl][5]))
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
                lat=row["lat"] or row["latitude"]
                long=row["long"]
                date=row["date"]
                print(lat,long,date)
                if lat == 'NA' and long == 'NA':
                    continue
                #else:
                    #long,lat=trans4mEPSG("EPSG:4326","EPSG:3035",float(long),float(lat))
                if float(lat) > self.minlat and float(lat) < self.maxlat and float(long) > self.minlong and float(long) < self.maxlong:
                    sampleinfo=(row["lat"],row["long"], row["date"])
                    filtered_data[row["sampleId"]] = sampleinfo
                    #filtered_data.append(sampleinfo)
        self.samples=filtered_data



class MemoryLogHandler(logging.Handler):
    """
    Custom logging handler that stores log messages in memory.
    """
    def __init__(self):
        super().__init__()
        self.logs = []  # List to store log records
    def emit(self, record):
        # Format the log message and add it to the list
        log_entry = self.format(record)
        self.logs.append(log_entry)

# Create a log object that uses the custom handler
class LogObject:
    def __init__(self):
        self.logger = logging.getLogger("MemoryLogger")
        self.logger.setLevel(logging.DEBUG)  # Set logging level
        # Add the custom handler
        self.memory_handler = MemoryLogHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        self.memory_handler.setFormatter(formatter)
        self.logger.addHandler(self.memory_handler)
    def get_logs(self):
        return self.memory_handler.logs  # Return the stored log messages
