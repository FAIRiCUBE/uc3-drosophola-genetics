import s3fs
import xarray as xr
from configparser import ConfigParser
import os

#Get UC1-Vienna100m data from Common S3 bucket (EOXHUB)

config_file="/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/S3/config.ini"
local_dir="/YOUR/PATH/S3/Vienna100m"

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


#def access_data_s3_bucket(filename: str) -> xr.DataArray:
#    """ Access a data cube from s3 bucket using xarray. Tries to use credentials stored in environment variables or in local config file """
#    if (os.environ.get('S3_FAIRICUBE_STORAGE_BUCKET') is None):
#        print('environment variables not set, trying with local config file')
#        try:
#            s3_config = config('/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/S3/config.ini')
#        except KeyError as er:
#            print(er)
#            print('Config file not found or malformed')
#    else:
#        s3_config = {
#        's3_fairicube_storage_bucket': os.environ.get('S3_FAIRICUBE_STORAGE_BUCKET'),
#        's3_fairicube_storage_key': os.environ.get('S3_FAIRICUBE_STORAGE_KEY'),
#        's3_fairicube_storage_secret': os.environ.get('S3_FAIRICUBE_STORAGE_SECRET')}
#    s3fs_FS = s3fs.S3FileSystem(
#        key=s3_config['s3_fairicube_storage_key'],
#        secret=s3_config['s3_fairicube_storage_secret'],
#    )
#    s3map = s3fs.S3Map(root=f's3:///{s3_config["s3_fairicube_storage_bucket"]}/{filename}', s3=s3fs_FS)
#    #ds = xr.open(store=s3map) # change this to the appropriate open_* function
#    ds = xr.open_zarr(fsspec.get_mapper(f's3:///{s3_config["s3_fairicube_storage_bucket"]}/{filename}', anon=True), consolidated=True)
#    return ds


###try again
def list_files_in_bucket(config_file: str):
    s3_config = config(config_file)
    
    # Extract S3 configuration
    aws_access_key_id = s3_config['s3_fairicube_storage_key']
    aws_secret_access_key = s3_config['s3_fairicube_storage_secret']
    region = s3_config['region']
    bucket_name = s3_config['s3_fairicube_storage_bucket']
    
    # Create S3 file system
    s3 = s3fs.S3FileSystem(
        key=aws_access_key_id,
        secret=aws_secret_access_key,
        client_kwargs={'region_name': region}
    )
    
    # List all items in the specified S3 bucket and directory
    #items = s3.ls(f's3://{bucket_name}/vienna_data/100m/', detail=True)
    items = s3.ls(f's3://{bucket_name}/vienna_data/100m/', detail=True)
    # Filter the items to include only files (not directories or the path itself)
    file_list = [item['Key'] for item in items if item['type'] == 'file' and not item['Key'].endswith('/')]
    dir_list = [item['Key'] for item in items if item['type'] == 'directory']
    #print("Files in bucket root:",file_list)
    #print("Directories in bucket root:", dir_list)
    return s3, file_list, dir_list  

#s3,files = list_files_in_bucket('/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/S3/config.ini')

def download_files_from_bucket(config_file: str, local_dir: str):
    s3, file_list, dir_list = list_files_in_bucket(config_file)
    # Ensure local directory exists
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)
    for file_path in file_list:
        # Remove the 's3://bucket_name/' prefix to get the file name
        file_name = file_path.split('/')[-1]
        print(file_path)
        local_file_path = os.path.join(local_dir, file_name)
        print(file_name)
        # Download the file
        s3.get(f's3://{file_path}', local_file_path)
        print(f"Downloaded {file_name} to {local_file_path}")

if __name__ == "__main__":
    download_files_from_bucket(config_file, local_dir)


### S4E Code
#if __name__ == "__main__":
#   filename = 'fairicube/vienna_data/100m/RR_2020-01-01_2024-05-01_epsg31256_100m_regridded.zarr'
#   ds = access_data_s3_bucket(filename)
#   print(ds)