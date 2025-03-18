import s3fs
import xarray as xr
from configparser import ConfigParser
import os

def config(filename: str, section: str='s3') -> dict:
    """ Helper function to parse .ini configuration files"""
    # create a parser
    parser = ConfigParser()
    # read config file
    parser.read(filename)
    # get section, default to postgresql
    section_dict = {}
    if parser.has_section(section):
        params = parser.items(section)
        for param in params:
            section_dict[param[0]] = param[1]
    else:
        raise KeyError(
            f'Section {section} not found in {filename}. Did you specifiy the right path?')
    return section_dict


def access_data_s3_bucket(filename: str) -> xr.DataArray:
    """ Access a data cube from s3 bucket using xarray. Tries to use credentials stored in environment variables or in local config file """
    if (os.environ.get('S3_FAIRICUBE_STORAGE_BUCKET') is None):
        print('environment variables not set, trying with local config file')
        try:
            s3_config = config('path/to/config.ini')
        except KeyError as er:
            print(er)
            print('Config file not found or malformed')
    else:
        s3_config = {
        's3_fairicube_storage_bucket': os.environ.get('S3_FAIRICUBE_STORAGE_BUCKET'),
        's3_fairicube_storage_key': os.environ.get('S3_FAIRICUBE_STORAGE_KEY'),
        's3_fairicube_storage_secret': os.environ.get('S3_FAIRICUBE_STORAGE_SECRET')}

    s3fs_FS = s3fs.S3FileSystem(
        key=s3_config['s3_fairicube_storage_key'],
        secret=s3_config['s3_fairicube_storage_secret'],
    )

    s3map = s3fs.S3Map(root=f's3:///{s3_config["s3_fairicube_storage_bucket"]}/{filename}', s3=s3fs_FS)
    ds = xr.open_zarr(store=s3map) # change this to the appropriate open_* function
    return ds

if __name__ == "__main__":
    filename = 'vienna_data/100m/TX_yearly_avg_2020-2024_epsg31256_100m_regridded.zarr'
    ds = access_data_s3_bucket(filename)
    print(ds)