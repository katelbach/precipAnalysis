import requests
import numpy as np
import pandas as pd
import datetime as dt
from html.parser import HTMLParser
from io import StringIO


class WyomingHTMLParser(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.data = None
        self.recording = 0

    def handle_starttag(self, tag, attributes):
        if tag != 'pre':
            return
        else:
            self.recording = 1

    def handle_endtag(self, tag):
        if tag == 'pre':
            self.recording = 0

    def handle_data(self, data):
        if self.recording:
            self.data = data


def get_radiosonde_data_from_html(text):

    parser = WyomingHTMLParser()
    parser.feed(text)
    parser.close()

    data = pd.read_table(
        StringIO(parser.data), skiprows=lambda x: x in [0, 1, 3, 4], sep=r"\s+")

    return data


def construct_url(time):

    url = f"http://weather.uwyo.edu/cgi-bin/bufrraob.py?datetime=" \
          f"{time:%Y-%m-%d %H:%M:%S}&id=11035&type=TEXT:LIST"

    return url


def get_profile(time):
    if time.minute != 0:
        raise ValueError('Radiosonde not available at minute: ', time.minute)

    if time.hour not in [0, 6, 12, 18]:
        raise ValueError('Radiosonde not available at hour: ', time.hour)

    url = construct_url(time)
    r = requests.get(url, allow_redirects=True)

    html = r.text

    data = get_radiosonde_data_from_html(html)

    return data


if __name__ == '__main__':
    time = dt.datetime(2024, 8, 17, 12, 0)
    data = get_profile(time)
