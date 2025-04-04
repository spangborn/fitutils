import argparse
from fitparse import FitFile
from geopy.distance import geodesic
from datetime import timedelta
import folium
import os


def extract_fit_data(fit_path):
    fitfile = FitFile(fit_path)
    gps_points = []
    power_data = []
    hr_data = []

    for record in fitfile.get_messages('record'):
        fields = {field.name: field.value for field in record}
        timestamp = fields.get('timestamp')

        if 'position_lat' in fields and 'position_long' in fields and timestamp:
            lat = fields['position_lat'] * (180 / 2**31)
            lon = fields['position_long'] * (180 / 2**31)
            gps_points.append((timestamp, lat, lon))

        if 'power' in fields and timestamp:
            power_data.append((timestamp, fields['power']))

        if 'heart_rate' in fields and timestamp:
            hr_data.append((timestamp, fields['heart_rate']))

    return gps_points, power_data, hr_data


def sync_using_multidata(base_power, compare_power, base_hr, compare_hr, base_gps, compare_gps):
    # Fallback if both power and HR data are missing, use GPS synchronization
    if not base_power and not base_hr and not compare_power and not compare_hr:
        print("Power and Heart Rate data are missing. Falling back to GPS synchronization.")
        return sync_using_gps(base_gps, compare_gps)  # Fall back to GPS sync if power and HR are missing

    min_offset = -15
    max_offset = 15
    best_offset = 0
    best_error = float('inf')

    base_power_dict = {t: p for t, p in base_power} if base_power else {}
    base_hr_dict = {t: hr for t, hr in base_hr} if base_hr else {}

    for offset in range(min_offset, max_offset + 1):
        error = 0
        matched = 0

        if compare_power:
            for t, p in compare_power:
                t_adj = t - timedelta(seconds=offset)
                # Only perform the calculation if a matching value is found in base_power_dict
                if t_adj in base_power_dict and base_power_dict[t_adj] is not None:
                    error += (base_power_dict[t_adj] - p) ** 2
                    matched += 1

        if compare_hr:
            for t, hr in compare_hr:
                t_adj = t - timedelta(seconds=offset)
                # Only perform the calculation if a matching value is found in base_hr_dict
                if t_adj in base_hr_dict and base_hr_dict[t_adj] is not None:
                    error += (base_hr_dict[t_adj] - hr) ** 2
                    matched += 1

        if matched > 0 and error < best_error:
            best_error = error
            best_offset = offset

    return best_offset


def sync_using_gps(base_gps, compare_gps):
    # Fallback sync method using only GPS data
    min_offset = -15
    max_offset = 15
    best_offset = 0
    best_error = float('inf')

    base_gps_dict = {t: (lat, lon) for t, lat, lon in base_gps}

    for offset in range(min_offset, max_offset + 1):
        error = 0
        matched = 0

        for t, lat, lon in compare_gps:
            t_adj = t - timedelta(seconds=offset)
            if t_adj in base_gps_dict:
                base_lat, base_lon = base_gps_dict[t_adj]
                error += geodesic((base_lat, base_lon), (lat, lon)).meters ** 2
                matched += 1

        if matched > 0 and error < best_error:
            best_error = error
            best_offset = offset

    return best_offset


def sync_tracks(base_track, compare_track, offset_seconds):
    shifted_compare = [(t - timedelta(seconds=offset_seconds), lat, lon) for t, lat, lon in compare_track]

    synced_base = []
    synced_compare = []

    base_dict = {t: (lat, lon) for t, lat, lon in base_track}
    for t, lat, lon in shifted_compare:
        closest_time = min(base_dict.keys(), key=lambda bt: abs((bt - t).total_seconds()))
        if abs((closest_time - t).total_seconds()) <= 30:
            synced_base.append((base_dict[closest_time][0], base_dict[closest_time][1]))
            synced_compare.append((lat, lon))

    return synced_base, synced_compare


def compare_tracks(base_track, compare_track):
    total_points = min(len(base_track), len(compare_track))
    if total_points == 0:
        return 0.0, 0.0, 0.0, 0, 0

    distances = []
    for base_point, compare_point in zip(base_track, compare_track):
        distances.append(geodesic(base_point, compare_point).meters)

    avg = sum(distances) / len(distances)
    min_dev = min(distances)
    max_dev = max(distances)
    min_index = distances.index(min_dev)
    max_index = distances.index(max_dev)

    return avg, min_dev, max_dev, min_index, max_index


def deviation_color(dist):
    if dist < 3:
        return "green"
    elif dist < 7:
        return "orange"
    else:
        return "red"


def large_emoji_div_icon(emoji):
    return folium.DivIcon(
        html=f"""<div style="font-size: 36px; transform: translate(-50%, -50%);">{emoji}</div>"""
    )


def plot_tracks(base_track, compare_track, timestamps, sync_index, min_index, max_index, title, output_name):
    sync_point = base_track[sync_index] if 0 <= sync_index < len(base_track) else None
    min_point = base_track[min_index] if 0 <= min_index < len(base_track) else None
    max_point = base_track[max_index] if 0 <= max_index < len(base_track) else None

    center_lat = (base_track[0][0] + base_track[-1][0]) / 2
    center_lon = (base_track[0][1] + base_track[-1][1]) / 2
    m = folium.Map(location=[center_lat, center_lon], zoom_start=16, tiles='OpenStreetMap')

    folium.PolyLine(base_track, color="blue", weight=4, opacity=0.7, tooltip="Base Track").add_to(m)
    folium.PolyLine(compare_track, color="red", weight=4, opacity=0.7, tooltip="Compare Track").add_to(m)

    last_time = None
    for t, base_point, compare_point in zip(timestamps, base_track, compare_track):
        if last_time is None or (t - last_time).total_seconds() >= 5:
            dist = geodesic(base_point, compare_point).meters
            color = deviation_color(dist)
            folium.CircleMarker(
                location=compare_point,
                radius=min(dist * 0.5, 12),
                color=color,
                fill=True,
                fill_opacity=0.5,
                tooltip=f"{dist:.1f} m @ {t.strftime('%H:%M:%S')}"
            ).add_to(m)
            last_time = t

    if sync_point:
        folium.Marker(sync_point, icon=large_emoji_div_icon("⭐"), tooltip="Sync Point").add_to(m)

    if min_point:
        folium.Marker(min_point, icon=large_emoji_div_icon("😊"), tooltip="Min Deviation").add_to(m)

    if max_point:
        folium.Marker(max_point, icon=large_emoji_div_icon("😭"), tooltip="Max Deviation").add_to(m)

    # Add legend
    legend_html = """
    <div style="position: fixed; bottom: 50px; left: 50px; width: 250px; height: 150px; background-color: white; 
    opacity: 0.7; z-index: 9999; border-radius: 10px; padding: 10px;">
        <h4>Deviation Legend</h4>
        <i style="background: green; width: 20px; height: 20px; display: inline-block; border-radius: 50%;"></i> Low Deviation<br>
        <i style="background: orange; width: 20px; height: 20px; display: inline-block; border-radius: 50%;"></i> Medium Deviation<br>
        <i style="background: red; width: 20px; height: 20px; display: inline-block; border-radius: 50%;"></i> High Deviation<br>
        <i style="font-size: 36px;">⭐</i> Sync Point<br>
        <i style="font-size: 36px;">😊</i> Min Deviation<br>
        <i style="font-size: 36px;">😭</i> Max Deviation<br>
    </div>
    """
    m.get_root().html.add_child(folium.Element(legend_html))

    output_file = f"{output_name.replace(' ', '_')}.html"
    m.save(output_file)
    print(f"Map saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Compare GPS waypoints from FIT files using power and heart rate data for sync.")
    parser.add_argument("files", nargs="+", help="FIT files to compare")
    args = parser.parse_args()

    if len(args.files) < 2:
        print("Please provide at least two FIT files to compare.")
        return

    all_data = {filename: extract_fit_data(filename) for filename in args.files}

    base_file = args.files[0]
    base_track, base_power, base_hr = all_data[base_file]

    print(f"\nBase file: {base_file}\n")

    for filename in args.files[1:]:
        other_track, other_power, other_hr = all_data[filename]

        offset = sync_using_multidata(base_power, other_power, base_hr, other_hr, base_track, other_track)
        synced_base, synced_other = sync_tracks(base_track, other_track, offset)

        if not synced_base or not synced_other:
            print(f"Skipping {filename} due to insufficient synced points.")
            continue

        avg_dev, min_dev, max_dev, min_idx, max_idx = compare_tracks(synced_base, synced_other)
        times = [t for t, _, _ in base_track]
        sync_index = len(synced_base) // 2

        print(f"Compared with {filename}:")
        print(f" - Time offset (via multi-data sync): {offset} seconds")
        print(f" - Synced GPS points: {len(synced_base)}")
        print(f" - Avg GPS deviation: {avg_dev:.2f} meters")
        print(f" - Min GPS deviation: {min_dev:.2f} meters at {times[min_idx]}")
        print(f" - Max GPS deviation: {max_dev:.2f} meters at {times[max_idx]}")

        plot_tracks(synced_base, synced_other, times, sync_index, min_idx, max_idx,
                    f"{os.path.basename(base_file)} vs {os.path.basename(filename)}",
                    f"Comparison_{os.path.basename(base_file)}_vs_{os.path.basename(filename)}")


if __name__ == "__main__":
    main()
