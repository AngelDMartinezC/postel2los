from sunpy.map import Map, make_fitswcs_header
from sunpy.util.exceptions import warn_metadata
from drms import Client

from astropy.coordinates import SkyCoord
import astropy.units as u

from datetime import datetime
import sys


def load_data(name):
    '''
    Read data into a sunpy.map.Map objecto.

    Parameters:
    ----------
    name: str
        Name of the frame to be reprojected

    Returns:
    -------
    sunpy.map.Map
        Sunpy object
    '''

    map_postel = Map(name)  # Load data

    # Check if T_REC keyword is defined. This keyword is used to extract the
    # coordinate information from the LOS header at that time, which will be
    # used as the template for the reprojection.
    try:
        map_postel.meta['T_REC']
    except KeyError:
        print('KeyError: \'T_REC\' keyword not found. Provide a date for '
              'T_REC in the header with the format YYYY.mm.dd_HH:MM:SS. '
              'Exiting.')
        sys.exit(1)

    # Check if RSUN_OBS is in the header. If not present, Sunpy will assign a
    # standard value for the solar radius. This might impair on the
    # reprojection.
    if 'RSUN_OBS' not in map_postel.meta:
        warn_metadata('Missing RSUN_OBS from metadata. The projection will '
                      'not be accurate unless the keyword RSUN_OBS is '
                      'specified.')

    return map_postel


def get_hmi_keywords(t_rec, obs='v'):
    '''
    Function get_hmi_keywords. Given the string date t_rec, in the form
    YYYY.mm.dd_HH:MM:SS, it extracts keywords to use them as template in the
    helioprojective header.

    Parameters:
    ----------
    t_rec : str
        Timestamp in the format YYYY.mm.dd_HH:MM:SS

    obs : str, optional
        Observing instrument identifier, default is 'v' for HMI. Other options
        are 'm' for LOS magnetic field, or 'ic' for continuum intensity.

    Returns:
    -------
    drms.QueryResponse or None
        Query response containing metadata keywords for reprojection.
        Returns None if there's an error during query execution.
    '''

    # Convert the T_REC time into a proper datetime class
    trecstep = 45  # TRECSTEP for HMI
    year = t_rec[0:4]
    month = t_rec[5:7]
    day = t_rec[8:10]
    hour = t_rec[11:13]
    minute = t_rec[14:16]
    second = t_rec[17:19]
    t_rec = f'{year}.{month}.{day}_{hour}:{minute}:{second}'
    time_0 = datetime.strptime(t_rec, r'%Y.%m.%d_%H:%M:%S')

    # Selects initial and final time
    qstr = f'hmi.{obs}_{int(trecstep)}s[{time_0}_TAI]'

    # Obtain coordinate metadata information
    client = Client(verbose=True)
    query = client.query(qstr, key='T_REC, CRPIX1, CRPIX2, CDELT1, CDELT2,'
                         'CRVAL1, CRVAL2, CROTA2, CRLN_OBS, CRLT_OBS,'
                         'TELESCOP, DSUN_OBS, RSUN_OBS, DATE-OBS,'
                         'CAR_ROT', convert_numeric=False)
    return query


def update_map_header(out_map, keywords, out_header):
    '''
    Update the header of out_map with keywords from the LOS template.

    Parameters:
    ----------
    out_map : sunpy.map.Map
        Reprojected solar map header to update.

    keywords : drms.QueryResponse
        Query response containing metadata keywords.

    out_header : dict
        Header dictionary containing updated metadata.
    '''

    out_map.meta['WCSNAME'] = 'Helioprojective-cartesian'

    # Update header with relevant keywords
    for key in ['CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2', 'CRVAL1', 'CRVAL2',
                'CROTA2', 'CRLN_OBS', 'CRLT_OBS', 'CAR_ROT', 'DSUN_OBS',
                'RSUN_OBS']:
        if key in keywords:
            out_map.meta[key] = float(keywords[key].values[0])

    out_map.meta['CUNIT1'] = out_header['CUNIT1']
    out_map.meta['CUNIT2'] = out_header['CUNIT2']
    out_map.meta['DATE-OBS'] = keywords['DATE-OBS'].values[0]

    # Keyword comments
    out_map.meta['keycomments'] = {
        'CDELT1': '[arcsec/pixel] image scale in the x direction',
        'CDELT2': '[arcsec/pixel] image scale in the y direction',
        'CRPIX1': '[pixel] CRPIX1: location of the Sun center',
        'CRPIX2': '[pixel] CRPIX2: location of the Sun center',
        'CRVAL1': '[arcsec] CRVAL1: x origin',
        'CRVAL2': '[arcsec] CRVAL2: y origin',
        'CRLN_OBS': '[deg] Carrington longitude of the observer',
        'CRLT_OBS': '[deg] Carrington latitude of the observer',
        'CROTA2': '[deg]',
        'CAR_ROT': 'Carrington rotation number of CRLN_OB',
        'DSUN_OBS': '[m] Distance from SDO to Sun center',
        'RSUN_OBS': '[arcsec] angular radius of Sun',
        'DATE-OBS': r'[ISO] Observation date {DATE__OBS}',
    }

    return None


# -----------------------------------------------------------------------


def main(name='INT.0000.fits', output='LOS_projection.fits'):
    '''
    Reproject a Postel-projected map into a LOS one.

    Parameters:
    ----------
    name: str
        Name of the Posel-projected map

    output: str
        Name of the reprojected LOS map
    '''
    # input_name = 'INT.0000.fits'
    # output_name = 'LOS_projection.fits'

    # Load data
    map_postel = load_data(name=name)

    # Shift the center to align with the center of the reprojected map
    map_postel.meta['CRPIX1'] += 1.5
    map_postel.meta['CRPIX2'] += 1.5

    # Get keywords from LOS map
    t_rec = map_postel.meta['T_REC']
    keywords = get_hmi_keywords(t_rec)
    naxis1 = 4096
    naxis2 = 4096
    crpix1 = float(keywords['CRPIX1'].values[0])
    crpix2 = float(keywords['CRPIX2'].values[0])
    cdelt1 = float(keywords['CDELT1'].values[0])
    cdelt2 = float(keywords['CDELT2'].values[0])
    crota2 = float(keywords['CROTA2'].values[0])

    # Retrieve data of center of map, then convert it to helioprojective
    crln_obs = map_postel.meta['CRLN_OBS']
    crlt_obs = map_postel.meta['CRLT_OBS']
    origin_hgcar = SkyCoord(crln_obs*u.deg, crlt_obs*u.deg, unit=u.deg,
                            frame=map_postel.coordinate_frame)
    origin = origin_hgcar.helioprojective

    # Create header of a helioprojective projection
    out_header = make_fitswcs_header(
        data=(naxis1, naxis2),  # size of output datacube
        scale=[cdelt1, cdelt2] * u.arcsec/u.pix,
        coordinate=origin,
        rotation_angle=crota2 * u.deg,
        projection_code="TAN",  # TAN: gnomonic projection
        reference_pixel=[crpix1 - 1, crpix2 - 1] * u.pix)

    # Reproject map from Postel into Helioprojective cartesian
    out_map = map_postel.reproject_to(out_header, algorithm='exact')

    # Upadate header of reprojected map
    update_map_header(out_map, keywords, out_header)

    # Save file
    out_map.save(output, overwrite=True)

    return None


main(name='INT.0000.fits', output='LOS_projection.fits')
