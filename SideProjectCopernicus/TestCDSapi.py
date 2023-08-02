### note, cdsapi requires key that is sored in an .cdsapirc file 
### pip install the follwing packages

import cdsapi
import datetime
import xarray as xr

def download_temperature_data(city, date):
    c = cdsapi.Client()
    # Define the data request parameters
    request = {
        'variable': '2m_temperature',
        'product_type': 'reanalysis',
        'year': date.year,
        'month': date.month,
        'day': date.day,
        'time': '00:00',
        'format': 'netcdf',
        'area': f"{city['lon']}/{city['lat']}/{city['lon']}/{city['lat']}",  # lon_min/lat_min/lon_max/lat_max
    }
    # Download the data
    filename = f"{city['name']}_{date.strftime('%Y%m%d')}_temperature.nc"
    c.retrieve('reanalysis-era5-single-levels', request, filename)

if __name__ == "__main__":
    city = {
        'name': 'Paris',
        'lon': 2.3522,
        'lat': 48.8566
    }
    date = datetime.datetime(2022, 8, 1)  # Replace with your desired date
    download_temperature_data(city, date)

####
# Replace with the actual path to your downloaded NetCDF file
nc_file_path = '/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/SideProjectCopernicus/Paris_20220801_temperature.nc'

# Read the NetCDF file using xarray
dataset = xr.open_dataset(nc_file_path)

# Access the temperature data variable
temperature_data = dataset['t2m']

# Now, temperature_data is an xarray DataArray containing the temperature values
print(temperature_data)

# Extract the temperature value as a float
temperature_value = float(temperature_data.values)

# Print the temperature value
print("Temperature (in Kelvin):", temperature_value, "K")

# Convert Kelvin to degrees Celsius
temperature_celsius = temperature_value - 273.15

# Print the result
print("Temperature in Celsius:", temperature_celsius, "°C")

###part two: comnined cities

import cdsapi
import datetime
import os

def download_temperature_data(cities, date, output_dir):
    c = cdsapi.Client()
    for city in cities:
        # Define the data request parameters
        request = {
            'variable': '2m_temperature',
            'product_type': 'reanalysis',
            'year': date.year,
            'month': date.month,
            'day': date.day,
            'time': '00:00',
            'format': 'netcdf',
            'area': f"{city['lon']}/{city['lat']}/{city['lon']}/{city['lat']}",  # lon_min/lat_min/lon_max/lat_max
        }
        # Download the data
        filename = f"{city['name']}_{date.strftime('%Y%m%d')}_temperature.nc"
        output_path = os.path.join(output_dir, filename)
        c.retrieve('reanalysis-era5-single-levels', request, output_path)

if __name__ == "__main__":
    cities = [
        {'name': 'Paris', 'lon': 2.3522, 'lat': 48.8566},
        {'name': 'London', 'lon': -0.1278, 'lat': 51.5074},
        {'name': 'New York', 'lon': -74.0060, 'lat': 40.7128}
    ]
    
    date = datetime.datetime(2022, 8, 1)  # Replace with your desired date
    output_directory = '/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/SideProjectCopernicus/'  # Replace with the directory where you want to save the data
    download_temperature_data(cities, date, output_directory)

import xarray as xr
import glob

def merge_temperature_data(input_dir, output_file):
    # Get a list of all NetCDF files in the input directory
    nc_files = glob.glob(os.path.join(input_dir, "*.nc"))
    # Open all NetCDF files and combine them into a single dataset
    datasets = [xr.open_dataset(file) for file in nc_files]
    combined_dataset = xr.concat(datasets, dim='city')
    # Save the combined dataset to a new NetCDF file
    combined_dataset.to_netcdf(output_file)
    # Close the individual datasets
    for dataset in datasets:
        dataset.close()

if __name__ == "__main__":
    input_directory = '/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/SideProjectCopernicus/'  # Replace with the directory containing the downloaded data
    output_file = '/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/SideProjectCopernicus//combined_temperature.nc'  # Replace with the desired output file path
    merge_temperature_data(input_directory, output_file)

# Read the NetCDF file using xarray
dataset = xr.open_dataset("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/SideProjectCopernicus//combined_temperature.nc")

# Access the temperature data variable
temperature_data = dataset['t2m']

# Now, temperature_data is an xarray DataArray containing the temperature values
print(temperature_data)

# Extract the temperature value as a float
temperature_value = float(temperature_data.values)

# Print the temperature value
print("Temperature (in Kelvin):", temperature_value, "K")