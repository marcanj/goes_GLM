import boto3
import botocore
from tabulate import tabulate
import datetime
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt


def generate_prefix(date_time):
    year = date_time.year
    day_of_year = date_time.timetuple().tm_yday
    hour = date_time.hour
    
    prefix = f'GLM-L2-LCFA/{year}/{day_of_year:03d}/{hour:02d}/'
    return prefix

def list_s3_files(bucket_name, prefix):
    # Create a boto3 client for S3
    s3 = boto3.client('s3', config=botocore.config.Config(signature_version=botocore.UNSIGNED))

    # List objects in the specified bucket and prefix
    response = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)

    # Check if the bucket contains objects
    if 'Contents' in response:
        files = response['Contents']
        
        # Extract and format file information
        table_data = []
        for file in files:
            file_info = [
                file['Key'],
                file['Size'],
                file['LastModified']
            ]
            table_data.append(file_info)
        
        print("Table found!")
    else:
        print("Table not found..")
    return table_data

def search_table(table_data, target_time):
    # Convert the given date_time to a datetime object
    day = target_time.day
    month = target_time.month
    
    for row in table_data:
        file_key = row[0]
        
        # Extract start and end times from the file key
        start_str = file_key.split('_s')[1].split('_e')[0]
        end_str = file_key.split('_e')[1].split('_c')[0]
        
        # Convert start and end times to datetime objects
        start_str_hour = start_str[7:9]
        start_str_min = start_str[9:11]
        start_str_sec = start_str[11:13]
        
        end_str_hour = end_str[7:9]
        end_str_min = end_str[9:11]
        end_str_sec = end_str[11:13]
        
        start_time = datetime.datetime(target_time.year, month, day, int(start_str_hour), int(start_str_min), int(start_str_sec))
        end_time = datetime.datetime(target_time.year, month, day, int(end_str_hour), int(end_str_min), int(end_str_sec))
        
        # Check if the target_time falls within the start and end times
        if start_time <= target_time <= end_time:
            return file_key

    return None

def get_file_key(bucket_name, date_time):
    
    prefix = generate_prefix(date_time)
    table_data = list_s3_files(bucket_name, prefix)
    file_key = search_table(table_data, date_time)
    return file_key
        

def download_and_read_goes16_data(bucket_name, file_key):
    s3 = boto3.resource('s3', config=botocore.config.Config(signature_version=botocore.UNSIGNED))
    try:
        s3.Bucket(bucket_name).download_file(file_key, '/tmp/' + file_key.split('/')[-1])
        return Dataset('/tmp/' + file_key.split('/')[-1], 'r')
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == '404':
            print("The object does not exist.")
        else:
            raise

# Define the bucket name and prefix (directory)
bucket_name = 'noaa-goes16'

# Get timestamp in the format YYYY-MM-DD HH:MM:SS
timestamp = '2024-05-03 06:35:00.000'

date_time = datetime.datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S.%f")

file_key = get_file_key(bucket_name, date_time)
print('file_key:', file_key)


nc_data = download_and_read_goes16_data(bucket_name, file_key)

specific_lat = 31.87347
specific_lon = -96.6622
tolerance = 0.5


# Check if data is loaded successfully
if nc_data:
    # Extract optical energy data
    flash_energy = nc_data.variables['flash_energy'][:]
    lat = nc_data.variables['flash_lat'][:]
    lon = nc_data.variables['flash_lon'][:]
    flash_time_offsets = nc_data.variables['flash_time_offset_of_first_event'][:]
    
    # Extract start time from the dataset attributes
    start_time_str = nc_data.getncattr('time_coverage_start')
    start_time = datetime.datetime.strptime(start_time_str, "%Y-%m-%dT%H:%M:%S.%fZ")

    # Convert time offsets to human-readable format
    flash_times = [start_time + datetime.timedelta(seconds=float(offset)) for offset in flash_time_offsets]

    indices = [
        i for i in range(len(lat))
        if abs(lat[i] - specific_lat) <= tolerance and abs(lon[i] - specific_lon) <= tolerance
    ]

    if indices:
        filtered_times = [flash_times[i] for i in indices]
        filtered_energy = [flash_energy[i] for i in indices]

        # Plot the optical energy related to lightning events over time
        plt.figure(figsize=(12, 6))
        plt.plot(filtered_times, filtered_energy, 'o-', label=f"Lat: {specific_lat}, Lon: {specific_lon}")
        plt.xlabel('Time')
        plt.ylabel('Optical Energy (J)')
        plt.title('Optical Energy of Lightning Events over Time')
        plt.legend()
        plt.grid(True)
        plt.show()
        
    
    # Calculate total time span
    if flash_times:
        total_time_span = max(flash_times) - min(flash_times)
        print(f"Total time span of plotted data: {total_time_span}")

    # Plot the optical energy related to lightning events
    plt.figure(figsize=(12, 6))
    sc = plt.scatter(lon, lat, c=flash_energy, cmap='viridis', s=20, marker='o', edgecolor='k')
    plt.colorbar(sc, label='Optical Energy (J)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Optical Energy of Lightning Events from GOES-16 GLM')
    plt.show()

    # Close the netCDF dataset
    nc_data.close()
else:
    print("Failed to load data.")

