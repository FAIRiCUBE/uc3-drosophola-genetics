import s3fs
import xarray as xr
from configparser import ConfigParser
import os

#Get UC1-Vienna100m data from Common S3 bucket (EOXHUB)

config_file="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/S3/config.ini"
#config_file="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/S3/s3.en4"
#local_dir="/YOUR/PATH/S3/Vienna100m"
local_dir="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/S3/ViennaData"

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



###try again
def list_files_in_bucket(config_file: str):
    s3_config = config(config_file)
    # Extract S3 configuration
    aws_access_key_id = s3_config['s3_fairicube_storage_key']
    aws_secret_access_key = s3_config['s3_fairicube_storage_secret']
    region = s3_config['region']
    bucket_name=s3_config['s3_fairicube_storage_bucket']
    # Create S3 file system
    s3 = s3fs.S3FileSystem(
        key=aws_access_key_id,
        secret=aws_secret_access_key,
        client_kwargs={'region_name': region}
    )
    items = s3.ls(f's3://{bucket_name}/vienna_data/100m/', detail=True)
    # Filter the items to include only files (not directories or the path itself)
    file_list = [item['Key'] for item in items if item['type'] == 'file' and not item['Key'].endswith('/')]
    dir_list = [item['Key'] for item in items if item['type'] == 'directory']
    return s3, file_list, dir_list  

#s3,files = list_files_in_bucket('/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/S3/config.ini')

def download_files_from_bucket(config_file: str, local_dir: str):
    s3, file_list, dir_list = list_files_in_bucket(config_file)
    # Ensure local directory exists
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)
    for file_path in dir_list:
        if not os.path.exists(file_path) and not file_path.endswith('.zarr'):
            newdir=os.path.join(local_dir,file_path)
            print("Make dir:", newdir)
            os.makedirs(newdir)
        else:
           newdir=local_dir
        items = s3.ls(f's3://{file_path}', detail=True)
        file_list = [item['Key'] for item in items if item['type'] == 'file' and not item['Key'].endswith('/')]
        for file in file_list: 
        #print(items)
            print(file)
            file_name = (file.split('/')[-2])+"/"+(file.split('/')[-1])
            print(file_name)
            local_file_path = os.path.join(newdir, file_name)
            print(local_file_path)
            s3.get(f's3://{file}', local_file_path)
        #print(f"Downloaded {file_name} to {local_file_path}")

if __name__ == "__main__":
    download_files_from_bucket(config_file, local_dir)


