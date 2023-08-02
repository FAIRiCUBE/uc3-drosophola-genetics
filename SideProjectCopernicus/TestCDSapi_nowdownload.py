import cdsapi
import datetime
import io
import json

def retrieve_temperature_data(city, date):
    c = cdsapi.Client()
    request = {
        'variable': '2m_temperature',
        'product_type': 'reanalysis',
        'year': date.year,
        'month': date.month,
        'day': date.day,
        'time': '00:00',
        'format': 'json',
        'area': f"{city['lon']}/{city['lat']}/{city['lon']}/{city['lat']}",  # lon_min/lat_min/lon_max/lat_max
    }
    # Retrieve the data as JSON
    response = c.retrieve('reanalysis-era5-single-levels', request)
    # Write the JSON response to a file-like object (StringIO)
    # Extract the temperature values from the JSON response
    # Extract the temperature values from the JSON response
    temperature_values = response['locations'][0]['data'][0]['value']
    # Close the client
    c.close()
    return temperature_values


if __name__ == "__main__":
    cities = [
        {'name': 'Paris', 'lon': 2.3522, 'lat': 48.8566},
        {'name': 'London', 'lon': -0.1278, 'lat': 51.5074},
        {'name': 'New York', 'lon': -74.0060, 'lat': 40.7128}
    ]
    date = datetime.datetime(2022, 8, 1)  # Replace with your desired date
    for city in cities:
        temperature_values = retrieve_temperature_data(city, date)
        print(f"Temperature values for {city['name']}: {temperature_values}")

###
import cdsapi
c = cdsapi.Client()
# Define the data request parameters
request = {
    'variable': '2m_temperature',
    'product_type': 'reanalysis',
    'year': 2022,
    'month': 8,
    'day': 1,
    'time': '00:00',
    'format': 'json',
    'area': '2.3522/48.8566/2.3522/48.8566',  # lon_min/lat_min/lon_max/lat_max
}
# Download the data
#filename = 'Paris_20230801_temperature.nc'
response=c.retrieve('reanalysis-era5-single-levels', request)
if __name__ == "__main__":
    download_temperature_data()

####
import cdsapi

