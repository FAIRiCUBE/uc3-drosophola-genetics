import requests

def get_copernicus_data():
    # Replace with the actual Copernicus API endpoint URL
    endpoint_url = "https://odp.dataspace.copernicus.eu/odata/v1" #https://documentation.dataspace.copernicus.eu/APIs/On-Demand%20Production%20API.html#:~:text=The%20ODP%20API%20endpoint%20is%20https%3A%2F%2Fodp.dataspace.copernicus.eu%2Fodata%2Fv1.%20The,endpoint%20supports%20both%20HTTP%20and%20HTTPS%20protocols.
    # Replace with your Copernicus API key or authentication token
    #endpoint_url="https://cds.climate.copernicus.eu/api/v2"
    api_key = "234647:f411d358-fa31-4df3-9cd5-9f26d4b63951"
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
    }
    try:
        response = requests.get(endpoint_url, headers=headers)
        if response.status_code == 200:
            # Successful response
            data = response.json()
            print(data)
        else:
            # Handle other status codes if needed
            print(f"Error: {response.status_code}, {response.text}")
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    get_copernicus_data()

https://my.cmems-du.eu/thredds/dodsC/global-reanalysis-phy-001-026-grepv1-monthly