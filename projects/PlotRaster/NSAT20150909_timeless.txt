File NSAT20150909_timeless.nc (NC_FORMAT_NETCDF4):

     2 variables (excluding dimension variables):
        8 byte int time[]   (Contiguous storage)  
            standard_name: time
            long_name: Time
            axis: T
            stored_direction: increasing
            type: double
            units: hours since 1900-01-01
            calendar: proleptic_gregorian
        float Tair[lon,lat]   (Chunking: [720,360])  (Compression: shuffle,level 9)
            _FillValue: 1.00000002004088e+20
            long_name: Near-Surface Air Temperature
            standard_name: air_temperature
            units: K
            coordinates: time lon lat

     2 dimensions:
        lon  Size:720 
            _FillValue: NaN
            standard_name: longitude
            units: degrees_east
            axis: X
            long_name: Longitude
            type: double
            valid_max: 360
            valid_min: -180
            positive: up
        lat  Size:360 
            _FillValue: NaN
            standard_name: latitude
            units: degrees_north
            axis: Y
            long_name: Latitude
            positive: up

    8 global attributes:
        title: WATCH Forcing Data methodology applied to ERA5 data
        institution: Copernicus Climate Change Service
        contact: http://copernicus-support.ecmwf.int
        comment: Methodology implementation for ERA5 and dataset production by B-Open Solutions for the Copernicus Climate Change Service in the context of contract C3S_25c
        Conventions: CF-1.7
        summary: ERA5 data regridded to half degree regular lat-lon; Genuine land points from CRU grid and ERA5 land-sea mask only; Tair elevation & bias-corrected using CRU TS4.04 mean monthly temperature and mean diurnal temperature range
        reference: Cucchi et al., 2020, Earth Syst. Sci. Data, 12(3), 2097–2120, doi:10.5194/essd-12-2097-2020; Weedon et al., 2014, Water Resources Res., 50, 7505-7514, doi:10.1002/2014WR015638; Harris et al., 2020, Scientific Data, 7(1), doi:10.1038/s41597-020-0453-3
        licence: The dataset is distributed under the Licence to Use Copernicus Products. The corrections applied are based upon CRU TS4.04, distributed under the Open Database License (ODbL)
