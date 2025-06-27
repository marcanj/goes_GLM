import boto3
import botocore
import datetime
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.dates as mdates
import pandas as pd
import geopandas as gpd
import contextily as ctx
from shapely.geometry import box, Point
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
import matplotlib
matplotlib.use('TkAgg')

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
    # Determine satellite longitude (GOES-16 or GOES-18)
    sat_lon = -75.0 if 'goes16' in bucket_name.lower() else -137.0

    # Approximate Earth radius and pixel size
    R_earth_km = 6371.0
    glm_pixel_km = 10.0  # Fixed conservative estimate

    # Convert km to degrees latitude
    dlat = glm_pixel_km / 110.574  # degrees per km latitude

    # Convert km to degrees longitude (account for latitude)
    dlon = glm_pixel_km / (111.320 * np.cos(np.radians(lat)))

    return dlat, dlon



def extract_filtered_flash_data(nc_data, start_time, target_time, lat0, lon0, tol, time_win):
    lat = nc_data.variables['flash_lat'][:]
    lon = nc_data.variables['flash_lon'][:]
    energy = nc_data.variables['flash_energy'][:]
    offsets = nc_data.variables['flash_time_offset_of_first_event'][:]

    abs_times = np.array([start_time + datetime.timedelta(seconds=float(t)) for t in offsets])
    
    indices = [
        i for i in range(len(lat))
        if abs(lat[i] - lat0) <= tol and abs(lon[i] - lon0) <= tol and
        abs((abs_times[i] - target_time).total_seconds()) <= time_win
    ]

    data = {
        'time': [abs_times[i] for i in indices],
        'lat': [lat[i] for i in indices],
        'lon': [lon[i] for i in indices],
        'energy': [energy[i] for i in indices],
    }
    return pd.DataFrame(data)


def extract_filtered_data(nc_data, var_prefix, start_time, target_time, lat0, lon0, tol, time_win):
    lat = nc_data.variables[f'{var_prefix}_lat'][:]
    lon = nc_data.variables[f'{var_prefix}_lon'][:]
    energy = nc_data.variables[f'{var_prefix}_energy'][:]
    offsets = nc_data.variables[f'{var_prefix}_time_offset'][:]
    
    abs_times = np.array([start_time + datetime.timedelta(seconds=float(t)) for t in offsets])

    indices = [
        i for i in range(len(lat))
        if abs(lat[i] - lat0) <= tol and
           abs(lon[i] - lon0) <= tol and
           abs((abs_times[i] - target_time).total_seconds()) <= time_win
    ]

    data = {
        'time': [abs_times[i] for i in indices],
        'energy': [energy[i] for i in indices],
        'lat': [lat[i] for i in indices],
        'lon': [lon[i] for i in indices],
        'type': [var_prefix] * len(indices)
    }

    return pd.DataFrame(data)

# ------------------ Plotting ------------------

def plot_energy_series(event_df, group_df, flash_df, target_time):
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot each dataset with a unique marker and label
    if not event_df.empty:
        ax.scatter(event_df['time'], event_df['energy'], label='Events', color='red', s=20, alpha=0.6)
    if not group_df.empty:
        ax.scatter(group_df['time'], group_df['energy'], label='Groups', color='blue', s=30, alpha=0.6, marker='x')
    if not flash_df.empty:
        ax.scatter(flash_df['time'], flash_df['energy'], label='Flashes', color='green', s=50, alpha=0.6, marker='^')

    # Format plot
    ax.axvline(target_time, color='black', linestyle='--', label='Target Time')
    ax.set_xlabel("Time")
    ax.set_ylabel("Energy (Joules)")
    ax.set_title("GLM Energy Time Series")
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()

def plot_group_line(group_df, target_time):
    if group_df.empty:
        print("No group data available for line plot.")
        return

    # Sort by time for proper line plotting
    group_df = group_df.sort_values('time')

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(group_df['time'], group_df['energy'], color='blue', marker='o', linestyle='-')
    ax.axvline(target_time, color='black', linestyle='--', label='Target Time')

    # Labels and formatting
    ax.set_title("Group Energy vs Time")
    ax.set_xlabel("Time")
    ax.set_ylabel("Energy (Joules)")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()
    
def plot_combined_spatial_and_energy(event_df, lat0, lon0, tol, bucket_name, pause_time=0.04):
    if event_df.empty:
        print("No events found.")
        return

    event_df = event_df.copy()
    event_df.sort_values('time', inplace=True)
    unique_times = sorted(event_df['time'].unique())

    energy_min = event_df['energy'].min()
    energy_max = event_df['energy'].max()
    dlat, dlon = 0.08, 0.084
    cmap = plt.cm.inferno
    norm = plt.Normalize(vmin=energy_min, vmax=energy_max)

    def project_to_3857(lat, lon):
        point = gpd.GeoSeries([Point(lon, lat)], crs="EPSG:4326").to_crs(epsg=3857).geometry[0]
        return point.x, point.y

    # Precompute cumulative energy sums grouped by time for time series subplot
    energy_sum = event_df.groupby('time')['energy'].sum().sort_index()
    times_sorted = energy_sum.index.to_list()
    integrated_energy_max = energy_sum.max()

    fig, (ax_map, ax_ts) = plt.subplots(1, 2, figsize=(16, 8))

    # Initialize colorbar once for map subplot
    pc_dummy = PatchCollection([], cmap=cmap, norm=norm, edgecolor='none', alpha=0.5)
    cbar = fig.colorbar(pc_dummy, ax=ax_map, label='Event Energy (J)')
    
    x0, y0 = project_to_3857(lat0, lon0)
    buffer_m = 2 * 111000 * tol  # rough degrees to meters

    for current_time in unique_times:
        # --- Spatial subplot ---
        ax_map.clear()
        patches_list, colors = [], []

        events_at_time = event_df[event_df['time'] == current_time]

        for _, row in events_at_time.iterrows():
            lat, lon, e = row['lat'], row['lon'], row['energy']
            x, y = project_to_3857(lat, lon)
            dx, dy = project_to_3857(lat + dlat / 2, lon + dlon / 2)
            dx -= x
            dy -= y
            rect = Rectangle((x - dx, y - dy), 2 * dx, 2 * dy)
            patches_list.append(rect)
            colors.append(e)

        if patches_list:
            pc = PatchCollection(patches_list, cmap=cmap, norm=norm, edgecolor='none', alpha=0.5)
            pc.set_array(np.array(colors))
            ax_map.add_collection(pc)

        # Reference location
        ax_map.plot(x0, y0, marker='*', color='red', markersize=15, label='Reference Location')
        ax_map.set_xlim(x0 - buffer_m, x0 + buffer_m)
        ax_map.set_ylim(y0 - buffer_m, y0 + buffer_m)
        ctx.add_basemap(ax_map, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:3857")

        ax_map.set_title(f'GLM Energy on Map\nTime: {current_time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]}')
        ax_map.legend()

        # --- Time series subplot ---
        ax_ts.clear()

        # Select all times <= current_time for cumulative plot
        idx = [t <= current_time for t in times_sorted]
        ax_ts.set_ylim(0, integrated_energy_max*1.1)
        ax_ts.plot(energy_sum.index, energy_sum.values, '--', color='blue', label='Cumulative Integrated Energy')
        ax_ts.axvline(current_time, color='black', linestyle='--', label='Current Time')

        # Set fixed vertical scale
        ax_ts.set_xlabel('Time (UTC)')
        ax_ts.set_ylabel('Cumulative Integrated Event Energy (Joules)')
        ax_ts.set_title('Integrated Event Energy over Time')
        ax_ts.legend()
        ax_ts.grid(True)
        fig.autofmt_xdate()

        plt.tight_layout()
        plt.pause(pause_time)

    plt.show()
    
def plot_event_sum(event_df, target_time):
    if event_df.empty:
        print("No event data available for sum plot.")
        return

    # Sum energy for exactly matching times
    energy_sum = event_df.groupby('time')['energy'].sum().sort_index()

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(energy_sum.index, energy_sum.values, 'o', color='blue', label='Event Energy integrated')
    ax.axvline(target_time, color='black', linestyle='--', label='Target Time')

    ax.set_xlabel('Time (UTC)')
    ax.set_ylabel('Integrated Event Energy (Joules)')
    ax.set_title('Integrated Event Energy near Target')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()


def plot_spatial_progression_on_map(event_df, lat0, lon0, tol, bucket_name, pause_time=0.01):
    if event_df.empty:
        print("No events found.")
        return

    event_df = event_df.copy()
    event_df.sort_values('time', inplace=True)
    unique_times = sorted(event_df['time'].unique())

    energy_min = event_df['energy'].min()
    energy_max = event_df['energy'].max()
    dlat, dlon = 0.08, 0.084
    cmap = plt.cm.inferno
    norm = plt.Normalize(vmin=energy_min, vmax=energy_max)

    # Prepare transformation to Web Mercator
    def project_to_3857(lat, lon):
        point = gpd.GeoSeries([Point(lon, lat)], crs="EPSG:4326").to_crs(epsg=3857).geometry[0]
        return point.x, point.y

    fig, ax = plt.subplots(figsize=(10, 10))
    pc = PatchCollection([], cmap=cmap, norm=norm, edgecolor='none')
    fig.colorbar(pc, ax=ax, label='Event Energy (J)')  # colorbar once

    for current_time in unique_times:
        ax.clear()
        patches_list, colors = [], []

        events_at_time = event_df[event_df['time'] == current_time]

        for _, row in events_at_time.iterrows():
            lat, lon, e = row['lat'], row['lon'], row['energy']
            x, y = project_to_3857(lat, lon)
            dx, dy = project_to_3857(lat + dlat / 2, lon + dlon / 2)
            dx -= x
            dy -= y
            rect = patches.Rectangle((x - dx, y - dy), 2 * dx, 2 * dy)
            patches_list.append(rect)
            colors.append(e)

        if patches_list:
            pc = PatchCollection(patches_list, cmap=cmap, norm=norm, edgecolor='none', alpha=0.5)
            pc.set_array(np.array(colors))
            ax.add_collection(pc)

        # Reference location
        x0, y0 = project_to_3857(lat0, lon0)
        ax.plot(x0, y0, marker='*', color='red', markersize=15, label='Reference Location')

        # Set bounds based on projected region
        buffer_m = 2 * 111000 * tol  # rough degrees to meters
        ax.set_xlim(x0 - buffer_m, x0 + buffer_m)
        ax.set_ylim(y0 - buffer_m, y0 + buffer_m)

        # Add map
        ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:3857")

        ax.set_title(f'GLM Energy on Map\nTime: {current_time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]}')
        ax.legend()
        plt.tight_layout()
        plt.pause(pause_time)
    plt.show(block=True) 
    plt.close()


def plot_group_progression(group_df, lat0, lon0, tol, pause_time=0.1):
    if group_df.empty:
        print("No groups found.")
        return

    group_df = group_df.copy()
    group_df.sort_values('time', inplace=True)
    unique_times = sorted(group_df['time'].unique())

    energy_min = group_df['energy'].min()
    energy_max = group_df['energy'].max()
    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=energy_min, vmax=energy_max)

    fig, ax = plt.subplots(figsize=(10, 8))
    sc = ax.scatter([], [], c=[], cmap=cmap, norm=norm)
    cbar = fig.colorbar(sc, ax=ax, label='Group Energy (J)')

    for current_time in unique_times:
        ax.clear()

        groups_at_time = group_df[group_df['time'] == current_time]
        lats = groups_at_time['lat'].values
        lons = groups_at_time['lon'].values
        energies = groups_at_time['energy'].values

        scatter = ax.scatter(lons, lats, c=energies, cmap=cmap, norm=norm, s=50, edgecolor='black')

        # Reference marker
        ax.plot(lon0, lat0, marker='*', color='red', markersize=15, label='Reference Location')

        ax.set_xlim(lon0 - 2*tol, lon0 + 2*tol)
        ax.set_ylim(lat0 - 2*tol, lat0 + 2*tol)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_title(f'GLM Group Centroids\nTime: {current_time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]}')
        ax.grid(True)
        ax.legend()
        plt.tight_layout()
        plt.pause(pause_time)

    plt.close()

def plot_all_groups(group_df, lat0, lon0, tol):
    if group_df.empty:
        print("No groups found.")
        return

    group_df = group_df.copy()

    # Filter by spatial tolerance (if needed)
    group_df = group_df[
        (np.abs(group_df['lat'] - lat0) <= tol) &
        (np.abs(group_df['lon'] - lon0) <= tol)
    ]

    if group_df.empty:
        print("No groups within specified tolerance.")
        return

    energy_min = group_df['energy'].min()
    energy_max = group_df['energy'].max()
    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=energy_min, vmax=energy_max)

    fig, ax = plt.subplots(figsize=(10, 8))

    lats = group_df['lat'].values
    lons = group_df['lon'].values
    energies = group_df['energy'].values

    sc = ax.scatter(lons, lats, c=energies, cmap=cmap, norm=norm, s=50, edgecolor='black')
    fig.colorbar(sc, ax=ax, label='Group Energy (J)')

    ax.plot(lon0, lat0, marker='*', color='red', markersize=15, label='Reference Location')

    ax.set_xlim(lon0 - 2*tol, lon0 + 2*tol)
    ax.set_ylim(lat0 - 2*tol, lat0 + 2*tol)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('All GLM Group Centroids')
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()


from matplotlib.animation import FuncAnimation
def plot_and_save_video(event_df, lat0, lon0, tol, bucket_name, filename='glm_energy_animation.mp4', fps=8):
    if event_df.empty:
        print("No events found.")
        return

    event_df = event_df.copy()
    event_df.sort_values('time', inplace=True)
    unique_times = sorted(event_df['time'].unique())

    energy_min = event_df['energy'].min()
    energy_max = event_df['energy'].max()
    dlat, dlon = 0.08, 0.084
    cmap = plt.cm.inferno
    norm = plt.Normalize(vmin=energy_min, vmax=energy_max)

    def project_to_3857(lat, lon):
        point = gpd.GeoSeries([Point(lon, lat)], crs="EPSG:4326").to_crs(epsg=3857).geometry[0]
        return point.x, point.y

    energy_sum = event_df.groupby('time')['energy'].sum().sort_index()
    integrated_energy_max = energy_sum.max()

    fig, (ax_map, ax_ts) = plt.subplots(1, 2, figsize=(18, 10))

    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    pc_dummy = PatchCollection([], cmap=cmap, norm=norm, edgecolor='none', alpha=0.5)
    cbar = fig.colorbar(pc_dummy, ax=ax_map, label='Event Energy (J)')

    x0, y0 = project_to_3857(lat0, lon0)
    buffer_m = 2 * 111000 * tol

    # We'll keep references to artists to update them instead of clearing all the time
    patches_list = []
    pc = None

    ref_loc_plot, = ax_map.plot(x0, y0, marker='*', color='red', markersize=15, label='Reference Location')
    ax_map.set_xlim(x0 - buffer_m, x0 + buffer_m)
    ax_map.set_ylim(y0 - buffer_m, y0 + buffer_m)
    ctx.add_basemap(ax_map, source=ctx.providers.OpenStreetMap.Mapnik, crs="EPSG:3857")

    ax_map.set_title('')
    ax_map.legend()

    ax_ts.set_ylim(0, integrated_energy_max * 1.1)
    ax_ts.set_xlabel('Time (UTC)')
    ax_ts.set_ylabel('Cumulative Integrated Event Energy (Joules)')
    ax_ts.set_title('Integrated Event Energy over Time')
    ax_ts.grid(True)

    vertical_line = ax_ts.axvline(unique_times[0], color='black', linestyle='--', label='Current Time')
    time_series_line, = ax_ts.plot([], [], '-', color='blue', label='Cumulative Integrated Energy')

    ax_ts.legend()
    fig.autofmt_xdate()

    def update(frame_idx):
        current_time = unique_times[frame_idx]
        for coll in ax_map.collections:
            coll.remove()

        events_at_time = event_df[event_df['time'] == current_time]

        patches_list = []
        colors = []

        for _, row in events_at_time.iterrows():
            lat, lon, e = row['lat'], row['lon'], row['energy']
            x, y = project_to_3857(lat, lon)
            dx, dy = project_to_3857(lat + dlat / 2, lon + dlon / 2)
            dx -= x
            dy -= y
            rect = Rectangle((x - dx, y - dy), 2 * dx, 2 * dy)
            patches_list.append(rect)
            colors.append(e)

        if patches_list:
            pc = PatchCollection(patches_list, cmap=cmap, norm=norm, edgecolor='none', alpha=0.5)
            pc.set_array(np.array(colors))
            ax_map.add_collection(pc)

        ax_map.set_title(f'GLM Energy on Map\nTime: {current_time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]}')
        ax_map.legend()

        # Update time series subplot
        idx = [t <= current_time for t in energy_sum.index]
        
        times_to_plot = np.array(energy_sum.index)[idx]
        values_to_plot = np.array(energy_sum.values)[idx]
        time_series_line.set_data(times_to_plot, values_to_plot)
        vertical_line.set_xdata(current_time)

        ax_ts.set_xlim(times_to_plot[0], times_to_plot[-1])

        return ax_map.collections + [ref_loc_plot, time_series_line, vertical_line]

    ani = FuncAnimation(fig, update, frames=len(unique_times), blit=False, repeat=False)

    # Save as video mp4 using ffmpeg writer
    print(f"Saving animation to {filename} ...")
    ani.save(filename, writer='ffmpeg', fps=fps, dpi=150)
    print("Save complete.")

    plt.close(fig)

# ------------------ MAIN EXECUTION ------------------

if __name__ == "__main__":
    bucket_name = 'noaa-goes19'
    timestamp = '2025-05-10 22:41:55.630'
    date_time = datetime.datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S.%f")
    specific_lat, specific_lon = 29.71028, -82.30736
    tolerance = 0.2
    time_window_sec = 0.4

    file_key = get_file_key(bucket_name, date_time)
    print('file_key:', file_key)

    nc_data = download_and_read_goes_data(bucket_name, file_key)

    if nc_data:
        start_str = nc_data.getncattr('time_coverage_start')
        start_time = datetime.datetime.strptime(start_str, "%Y-%m-%dT%H:%M:%S.%fZ")
        target_time = date_time

        event = extract_filtered_data(nc_data, 'event', start_time, target_time, specific_lat, specific_lon, tolerance, time_window_sec)
        group = extract_filtered_data(nc_data, 'group', start_time, target_time, specific_lat, specific_lon, tolerance, time_window_sec)
        flash = extract_filtered_flash_data(nc_data, start_time, target_time, specific_lat, specific_lon, tolerance, time_window_sec)


        #plot_energy_series(event, group, flash, target_time)
        #plot_group_line(group, target_time)
        #plot_event_sum(event, target_time)
        #plot_spatial_progression_on_map(event, specific_lat, specific_lon, tolerance, bucket_name)
        #plot_combined_spatial_and_energy(event, specific_lat, specific_lon, tolerance, bucket_name)
        plot_and_save_video(event, specific_lat, specific_lon, tolerance, bucket_name, filename='glm_energy_animation.mp4', fps=8)

        #plot_group_progression(group, specific_lat, specific_lon, tolerance)
        #plot_all_groups(group, specific_lat, specific_lon, tolerance)
        nc_data.close()
    else:
        print("Failed to load data.")
