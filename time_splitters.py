"""

algorithm to convert from Julian Day to Gregorian date
Taken from:
Richards, E. G. (2013). Calendars. In S. E. Urban & P. K. Seidelmann, eds. Explanatory Supplement to the Astronomical
Almanac, 3rd ed. (pp. 585–624). Mill Valley, Calif.: University Science Books. pp617–9
"""

from math import floor
import numpy as np
from datetime import datetime, timedelta
from jdcal import gcal2jd

def JDayToGreg(J, ms):
    y = 4716
    j = 1401
    m = 2
    n = 12
    r = 4
    p = 1461
    v = 3
    u = 5
    s = 153
    w = 2
    B = 274277
    C = -38
    f = J + j + floor((floor((4*J+B)/146097)*3)/4) + C
    e = r*f+v
    g = floor((e%p)/r)
    h = u*g+w
    D = floor((h%s)/u) + 1
    M = (floor(h/s)+m)%n + 1
    Y = floor(e/p) - y + floor((n+m-M)/n)
    GDay = datetime(Y, M, D)
    delta = timedelta(milliseconds=int(ms))
    timestamp = GDay + delta
    return timestamp


def getLocalTS(ts, lon):
    difference_in_hours = timedelta(hours = (lon/360)*24)
    ts = ts + difference_in_hours
    if lon > 180:
        ts = ts - timedelta(hours = 24)
    tsf = ts.strftime("%Y_%m_%d_%H%M")
    return tsf

# get local time from GMT (UTC) time and longitude
# this is neither exact solar time nor geopolitically defined regional time, for the following reasons:
# it's based on the recorded GMT time, which is not solar
# packages like timezonefinder don't have timezones for the middle of the ocean, so can't be relied on for this
# ... so this is a compromise
def get_local_time(time, lon):
    difference_in_hours = get_local_timediff(lon)
    local_time = time + difference_in_hours
    return local_time

def get_local_timediff(lon):
    difference_in_hours = timedelta(hours = (lon/360)*24)
    if lon > 180:                           # if we're more than half way round the globe...
        difference_in_hours -= timedelta(hours = 24) # ...it's the day before
    return difference_in_hours

# takes a datetime.datetime object; returns julian day and ms since start of that julian day
def datetime_to_jday(datetime):
    year_greg = datetime.year
    month_greg = datetime.month
    day_greg = datetime.day
    jday = int(sum(gcal2jd(year_greg, month_greg, day_greg)) + 0.5)  # convert to Julian Day
    return jday

def datetime_to_ms(datetime):
    hours = datetime.hour
    minutes = datetime.minute
    seconds = datetime.second
    ms = (hours * 3600 + minutes * 60 + seconds) * 1000
    return ms

