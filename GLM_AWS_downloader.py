import boto3
import botocore
import datetime
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.dates as mdates
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

# ------------------ S3 Utilities ------------------

def generate_prefix(date_time):
    year = date_time.year
    day_of_year = date_time.timetuple().tm_yday
    hour = date_time.hour
    return f'GLM-L2-LCFA/{year}/{day_of_year:03d}/{hour:02d}/'

def list_s3_files(bucket_name, prefix):
    s3 = boto3.client('s3', config=botocore.config.Config(signature_version=botocore.UNSIGNED))
    response = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
    if 'Contents' in response:
        print("Table found!")
        return [[f['Key'], f['Size'], f['LastModified']] for f in response['Contents']]
    print("Table not found..")
    return []

def search_table(table_data, target_time):
    day, month = target_time.day, target_time.month
    for row in table_data:
        key = row[0]
        try:
            start_str = key.split('_s')[1].split('_e')[0]
            end_str = key.split('_e')[1].split('_c')[0]
            start = datetime.datetime(target_time.year, month, day,
                                      int(start_str[7:9]), int(start_str[9:11]), int(start_str[11:13]))
            end = datetime.datetime(target_time.year, month, day,
                                    int(end_str[7:9]), int(end_str[9:11]), int(end_str[11:13]))
            if start <= target_time <= end:
                return key
        except Exception:
            continue
    return None

def get_file_key(bucket_name, date_time):
    prefix = generate_prefix(date_time)
    table_data = list_s3_files(bucket_name, prefix)
    return search_table(table_data, date_time)

def download_and_read_goes_data(bucket_name, file_key):
    s3 = boto3.resource('s3', config=botocore.config.Config(signature_version=botocore.UNSIGNED))
    local_file = '/tmp/' + file_key.split('/')[-1]
    try:
        s3.Bucket(bucket_name).download_file(file_key, local_file)
        return Dataset(local_file, 'r')
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == '404':
            print("The object does not exist.")
        else:
            raise

# ------------------ Processing ------------------

def estimate_glm_pixel_deg(lat, lon, bucket_name):
    sat_lon = -75.0 if 'goes16' in bucket_name.lower() else -137.0
    def haversine(lat1, lon1, lat2, lon2):
        R = 6371
        phi1, phi2 = np.radians(lat1), np.radians(lat2)
        dphi = np.radians(lat2 - lat1)
        dlambda = np.radians(lon2 - lon1)
        a = np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2
        return 2 * R * np.arcsin(np.sqrt(a))
    distance_km = haversine(0.0, sat_lon, lat, lon)
    pixel_km = np.clip(8.0 + (distance_km / 5000) * (14.0 - 8.0), 8.0, 14.0)
    return pixel_km / 111.0, pixel_km / (111.0 * np.cos(np.radians(lat)))

def extract_filtered_data(nc_data, var_prefix, start_time, target_time, lat0, lon0, tol, time_win):
    lat = nc_data.variables[f'{var_prefix}_lat'][:]
    lon = nc_data.variables[f'{var_prefix}_lon'][:]
    energy = nc_data.variables[f'{var_prefix}_energy'][:]
    offsets = nc_data.variables[f'{var_prefix}_time_offset'][:]
    abs_times = np.array([start_time + datetime.timedelta(seconds=float(t)) for t in offsets])
    indices = [
        i for i in range(len(lat))
        if abs(lat[i] - lat0) <= tol and abs(lon[i] - lon0) <= tol and
        abs((abs_times[i] - target_time).total_seconds()) <= time_win
    ]
    return [abs_times[i] for i in indices], [energy[i] for i in indices]

# ------------------ Plotting ------------------

def plot_energy_series(event, group, flash, target_time):
    plt.figure(figsize=(12, 6))
    plt.plot(*event, 'o', label='Event Energy', color='blue')
    plt.plot(*group, 'x', label='Group Energy', color='green')
    plt.plot(*flash, '^', label='Flash Energy', color='red')
    plt.axvline(target_time, color='black', linestyle='--', label='Target Time')
    plt.xlabel('Time (UTC)')
    plt.ylabel('Optical Energy (Joules)')
    plt.title('GLM Optical Energy near Target Location and Time')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_group_line(group, target_time):
    if group[0] and group[1]:
        idx = np.argsort(group[0])
        sorted_times = [group[0][i] for i in idx]
        sorted_energies = [group[1][i] for i in idx]
        plt.figure(figsize=(12, 4))
        plt.plot(sorted_times, sorted_energies, 'o', color='green', label='Group Energy')
        plt.axvline(target_time, color='black', linestyle='--', label='Target Time')
        plt.xlabel('Time (UTC)')
        plt.ylabel('Group Optical Energy (Joules)')
        plt.title('Group Energy Over Time Near Target')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

def plot_event_sum(event, target_time):
    if event[0] and event[1]:
        energy_by_time = defaultdict(float)
        for t, e in zip(event[0], event[1]):
            energy_by_time[t] += e
        times = sorted(energy_by_time.keys())
        energies = [energy_by_time[t] for t in times]
        plt.figure(figsize=(12, 4))
        plt.plot(times, energies, 'o', color='blue', label='Event Energy integrated')
        plt.axvline(target_time, color='black', linestyle='--', label='Target Time')
        plt.xlabel('Time (UTC)')
        plt.ylabel('Integrated Event Energy (Joules)')
        plt.title('Integrated Event Energy near Target')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

def plot_spatial(nc_data, lat0, lon0, tol, bucket_name):
    lats = nc_data.variables['event_lat'][:]
    lons = nc_data.variables['event_lon'][:]
    energies = nc_data.variables['event_energy'][:]
    dlat, dlon = estimate_glm_pixel_deg(lat0, lon0, bucket_name)

    fig, ax = plt.subplots(figsize=(10, 8))
    patches_list, colors = [], []

    for lat, lon, e in zip(lats, lons, energies):
        if abs(lat - lat0) <= tol and abs(lon - lon0) <= tol:
            rect = patches.Rectangle((lon - dlon/2, lat - dlat/2), dlon, dlat)
            patches_list.append(rect)
            colors.append(e)

    if patches_list:
        pc = PatchCollection(patches_list, cmap='hot', edgecolor='none')
        pc.set_array(np.array(colors))
        ax.add_collection(pc)
        plt.colorbar(pc, ax=ax, label='Event Energy (J)')

        # 🔴 Add marker for reference location
        ax.plot(lon0, lat0, marker='*', color='red', markersize=15, label='Reference Location')

        ax.set_xlim(lon0 - 2*tol, lon0 + 2*tol)
        ax.set_ylim(lat0 - 2*tol, lat0 + 2*tol)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_title('GLM Event Energy Footprints')
        ax.grid(True)
        ax.legend()
        plt.tight_layout()
        plt.show()
    else:
        print("No events found for spatial plot.")


# ------------------ MAIN EXECUTION ------------------

if __name__ == "__main__":
    bucket_name = 'noaa-goes19'
    timestamp = '2025-05-10 22:41:55.630'
    date_time = datetime.datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S.%f")
    specific_lat, specific_lon = 29.71028, -82.30736
    tolerance = 0.2
    time_window_sec = 20

    file_key = get_file_key(bucket_name, date_time)
    print('file_key:', file_key)

    nc_data = download_and_read_goes_data(bucket_name, file_key)

    if nc_data:
        start_str = nc_data.getncattr('time_coverage_start')
        start_time = datetime.datetime.strptime(start_str, "%Y-%m-%dT%H:%M:%S.%fZ")
        target_time = date_time

        event = extract_filtered_data(nc_data, 'event', start_time, target_time, specific_lat, specific_lon, tolerance, time_window_sec)
        group = extract_filtered_data(nc_data, 'group', start_time, target_time, specific_lat, specific_lon, tolerance, time_window_sec)

        flash_lat = nc_data.variables['flash_lat'][:]
        flash_lon = nc_data.variables['flash_lon'][:]
        flash_energy = nc_data.variables['flash_energy'][:]
        flash_offsets = nc_data.variables['flash_time_offset_of_first_event'][:]
        flash_times = np.array([start_time + datetime.timedelta(seconds=float(t)) for t in flash_offsets])
        flash_indices = [
            i for i in range(len(flash_lat))
            if abs(flash_lat[i] - specific_lat) <= tolerance and
               abs(flash_lon[i] - specific_lon) <= tolerance and
               abs((flash_times[i] - target_time).total_seconds()) <= time_window_sec
        ]
        flash = ([flash_times[i] for i in flash_indices], [flash_energy[i] for i in flash_indices])

        plot_energy_series(event, group, flash, target_time)
        plot_group_line(group, target_time)
        plot_event_sum(event, target_time)
        plot_spatial(nc_data, specific_lat, specific_lon, tolerance, bucket_name)

        nc_data.close()
    else:
        print("Failed to load data.")
