#!/usr/bin/env python
'''
Python script to properly unstack a frame from a datacube (first argument).
The script computes the time of the frame (second argument) and extract the
observational Carrington Longitude and Latitude for that time. These
coordinates, as well as the time and several other keywords, are then put in
the unstacked frame.

    Usage: unstack_hmi.py input[FITS] frame_number[INT] suffix

'''
from datetime import datetime, timedelta
from astropy.io import fits
import drms
import sys
import os
from numpy import pi

# Check that the number of arguments is right
if len(sys.argv) != 4:
    sys.tracebacklimit = 0
    num = len(sys.argv)
    raise TypeError(f'unstack_hmi.py takes 3 arguments ({num-1} given)\n\n'
                    '  Usage: unstack_hmi.py input[FITS] frame_number[INT] '
                    'suffix\n')

# Read input variables
name = sys.argv[1]
frame = int(sys.argv[2])
prefix = sys.argv[3]

# Read datacube information
hdul = fits.open(name)
header = hdul[0].header
data = hdul[0].data

# Read keywords of interest
trecstep = header['DAXIS3']
t_rec = header['T_REC']
crpix1 = header['X0']
crpix2 = header['Y0']
crval1 = header['REF_L0']
crval2 = header['REF_B0']
daxis1 = header['DAXIS1']  # Custom keywords (units of solar radius)
daxis2 = header['DAXIS2']  # Custom keywords (units of solar radius)


def get_hmi_keywords(t_rec, obs='v'):
    '''
    Function get_hmi_keywords. Given a string date, t_rec, having the form
    YYYY.mm.dd_HH:MM:SS, reads useful keywords to be put in the unstacked
    header. This properly assings coordinates for unstacked frames from
    datacubes by reading keywords at the time of the unstacked frame.
    '''
    # Convert t_rec time into proper datetime class
    year = t_rec[0:4]
    month = t_rec[5:7]
    day = t_rec[8:10]
    hour = t_rec[11:13]
    minute = t_rec[14:16]
    second = t_rec[17:19]
    t_rec = f'{year}.{month}.{day}_{hour}:{minute}:{second}'
    print(f't_rec is: {t_rec}_TAI')
    time_0 = datetime.strptime(t_rec, r'%Y.%m.%d_%H:%M:%S')

    # Get TAI time of selected frame and convert it to string
    delta_time = timedelta(seconds=trecstep*frame)
    frame_date = time_0 + delta_time
    frame_date = datetime.strftime(frame_date, r'%Y.%m.%d_%H:%M:%S')

    # Selects initial and final time
    qstr = f'hmi.{obs}_{int(trecstep)}s[{frame_date}_TAI]'

    print('Reading metadata information from JSOC')
    # Obtain coordinate metadata information
    client = drms.Client(verbose=True)
    query = client.query(qstr, key='T_REC, CRPIX1, CRPIX2, CDELT1, CDELT2,'
                         'CRVAL1, CRVAL2, CROTA2, CRLN_OBS, CRLT_OBS,'
                         'TELESCOP, DSUN_OBS, RSUN_OBS, DATE-OBS,'
                         'T_REC_step, CAR_ROT', convert_numeric=False)
    return query


# Read info of keywords for the time specified by t_rec
query = get_hmi_keywords(t_rec)

# Unstack one frame and prepare it for output
unstacked = data[frame]
hdu_out = fits.PrimaryHDU(unstacked)
hdr_out = hdu_out.header

# Assing values to the keywords of the output header
hdr_out['WCSNAME'] = ('ARC', 'Postel projection')
hdr_out['CTYPE1'] = ('CRLN-ARC', 'x-axis tangent to constant longitude')
hdr_out['CTYPE2'] = ('CRLT-ARC', 'y-axis tangent to constant latitude')
hdr_out['CRPIX1'] = (crpix1, 'x-pix location of projection center')
hdr_out['CRPIX2'] = (crpix2, 'y-pix location of projection center')
hdr_out['CRVAL1'] = (crval1, '[deg] Longitude of projection at X0/CRPIX1')
hdr_out['CRVAL2'] = (crval2, '[deg] Latitude of projection at Y0/CRPIX2')
hdr_out['CDELT1'] = (180*daxis1/pi, '[deg] Arc len of x-pix separation')
hdr_out['CDELT2'] = (180*daxis2/pi, '[deg] Arc len of y-pix separation')
hdr_out['CUNIT1'] = ('deg', 'x-axis units are degrees of arc')
hdr_out['CUNIT2'] = ('deg', 'y-axis units are degrees of arc')
hdr_out['CRLN_OBS'] = (float(query['CRLN_OBS'].values[0]),
                       '[deg] Longitude at Solar Disk center')
hdr_out['CRLT_OBS'] = (float(query['CRLT_OBS'].values[0]),
                       '[deg] Latitude at Solar Disk center')
hdr_out['DSUN_OBS'] = (float(query['DSUN_OBS'].values[0]),
                       '[m] Distance from SDO to Sun center')
hdr_out['RSUN_OBS'] = (float(query['RSUN_OBS'].values[0]),
                       '[arcsec] angular radius of Sun')
hdr_out['TELESCOP'] = (query['TELESCOP'].values[0], 'Telescope')
hdr_out['TRECSTEP'] = (float(query['T_REC_step'].values[0]), '[seconds]')
hdr_out['T_REC'] = (query['T_REC'].values[0], '[TAI] Slot time')
hdr_out['DATE-OBS'] = (query['DATE-OBS'].values[0],
                       r'[ISO] Observation date {DATE__OBS}')

# save into file
name = f'{prefix}.{frame:04d}.fits'
if os.path.isfile(name):
    overwrite = input(f'File named {name} already exists.  Overwrite? [y/n] ')
    if (overwrite == 'y') or (overwrite == 'Y'):
        hdu_out.writeto(name, overwrite=True)
    else:
        print('File not overwitten')
else:
    hdu_out.writeto(name)

