import pandas as pd
import datetime as dt
import numpy as np
import json
import time
from pathlib import Path
from urllib.request import urlretrieve
import urllib.error

from pyextremes import EVA


def read_csv_data():
    columns = ["time", "rr"]

    df = pd.read_csv("data/PRECIP_1HR_VIENNA_1941_2024.csv", usecols=columns,
                     parse_dates=[0], index_col=0)

    df['2H_SUM'] = df['rr'].rolling('2h').sum()

    return df

def read_json_data(filepath):

    with open(filepath) as json_file:
        jsondata = json.load(json_file)['features']

    rr_data = dict()

    for station in jsondata:
        station_id = station['properties']['station']
        lat, lon = station['geometry']['coordinates']

        if station['properties']['parameters']['rr']['data'][0] is None:
            continue

        rr_sum = np.sum(station['properties']['parameters']['rr']['data'])
        rr_data[station_id] = {'lat': lat, 'lon': lon, 'rr': rr_sum}

    return rr_data


def get_json_data(time0, time1):

    station_ids = [105, 5925, 5802, 5935, 107, 5820, 106, 4115, 51, 222, 31,
                   5990]

    api_url = (f"https://dataset.api.hub.geosphere.at/v1/station/historical/"
               f"klima-v2-1h?parameters=RR"
               f"&station_ids={','.join([str(x) for x in station_ids])}"
               f"&start={(time0+dt.timedelta(hours=1)):%Y-%m-%dT%H:%M}"
               f"&end={time1:%Y-%m-%dT%H:%M}")

    filepath = Path("data", f"{time0:%Y%m%d}")
    filepath.mkdir(exist_ok=True)
    filepath = Path(filepath, "tawes.json")

    remaining_download_tries = 5
    while remaining_download_tries > 0:
        try:
            print("fetching TAWES data from GeoSphere datahub")
            urlretrieve(api_url, filepath)
            break
        except urllib.error.HTTPError:
            remaining_download_tries -= 1
            time.sleep(0.8)
            print("retrying...")
            continue

    if remaining_download_tries == 0:
        raise urllib.error.HTTPError
    return


def evt_analysis(df):

    # for the EVT analysis we exclude the event from 2024
    df_subset = df['1941-01-01': '2024-08-01']

    series = (
        df_subset['2H_SUM']
        .sort_index(ascending=True)
        .astype(float)
        .dropna()
    )

    model = EVA(series)
    model.get_extremes(method="BM", block_size="365.2425D")
    model.fit_model(distribution='genextreme')

    #model.plot_diagnostic(alpha=0.95)
    # model.plot_extremes()
    #plt.show()

    return model


