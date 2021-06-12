import geojson
import json
import numpy as np
import os
import pyproj
import rasterio
import shapely
import shapely.wkt
import subprocess

from osgeo import gdal
from rasterio.features import sieve, shapes
from scipy.ndimage import gaussian_filter
from shapely.ops import transform
from shapely.geometry import shape
from statsmodels.distributions.empirical_distribution import ECDF
        
def get_utm_zone(latitude, longitude):
    """
    compute utm zone and weather it is in North or South given by a lat/lon coordinate

    Arguments
      longitude : float
      latitude : float

    Returns
      utm_zone, is_north : list (or list like)
      utm zone number and N or S string

    """

    utm_zone = int(1 + (longitude + 180.0) / 6.0)

    is_north = 0
    if (latitude < 0.0):
        is_north = "S";
    else:
        is_north = "N";

    return utm_zone, is_north

def convert_wkt_to_geometry(geometry_wkt):
    """ 
    Convert wkt string to a shapely.geometry.polygon.Polygon object
    
    Arguments
      geometry_wkt : wkt string

    Returns
      geom: shapely geometry object

    """
    
    geom = shapely.wkt.loads(geometry_wkt)

    return geom


def compute_centroid_from_geometry(geometry_wkt):
    """
    compute centroid of a geometry; can be polygon, point

    Arguments
      geometry : str
      geojson geometry string

    Returns
      y, x: latitude and longitude of centroid

    """

    geometry = shapely.wkt.loads(geometry_wkt)
    x = geometry.centroid.x
    y = geometry.centroid.y

    return y, x

def convert_geojson_to_wkt(boundary):
    """
    Returns wkt geometry from geojson 

    Arguments:
    ----------
    geojson : json
        geojson 

    Returns:
    -------
    geometry wkt : wkt 
    """    
    
    s = json.dumps(boundary)
    # convert to geojson.geometry.Polygon
    geo = geojson.loads(s)
    # Feed to shape() to convert to shapely.geometry.polygon.Polygon
    geom = shape(geo)
    # get a WKTrepresentation
    return geom.wkt
    
def convert_geom_latlon_to_utm_na(geometry, utmzone):
    """
    reproject geometry in wgs84 lat/lon to utm projection for nothern hemisphere
    
    Arguments
    geometry: wkt string
        wkt format geometry string
    utmzone: string or int
        2 digit utmzone
    
    Return: None 
    """
    source_crs = 'epsg:4326'
    target_crs = 'epsg:326{}'.format(utmzone)

    project = partial(
        pyproj.transform,
        pyproj.Proj(init = source_crs), # source coordinate system
        pyproj.Proj(init = target_crs)) # destination coordinate system

    utm_geom = transform(project, geometry)  # apply projection
    return utm_geom

def convert_geom_latlon_to_utm_sa(geometry, utmzone):
    """
    reproject geometry in wgs84 lat/lon to utm projection for southern hemisphere
    
    Arguments
    geometry: wkt string
        wkt format geometry string
    utmzone: string or int
        2 digit utmzone
    
    Return: None 
    """

    source_crs = 'epsg:4326'
    target_crs = 'epsg:327{}'.format(utmzone)

    project = partial(
        pyproj.transform,
        pyproj.Proj(init = source_crs), # source coordinate system
        pyproj.Proj(init = target_crs)) # destination coordinate system

    utm_geom = transform(project, geometry)  # apply projection
    return utm_geom


def gdal_convert_raster_to_shp(input_raster, output_shp):
    """
    convert raster file to ESRI shapefile
    
    Arguments
    input_raster: string
        path to a input raster file
    
    Return: None 
    """
    cmd = "gdal_polygonize.py " + \
          input_raster + " " + \
          output_shp
    subprocess.check_call(cmd, shell=True)
    
def read_mask_image(raster_path, mask_feature):
    """
    clipping (masking) raster data to polygon geometry 
    
       
    Arguments
    raster_path: string
        path to the raster file
    
    mask_feature: shapely geometry object
        * be sure to be in list format

    Returns
    arr_ma: nd array
        clipped nd array
    arr_ma_tranform: rasterio geotransform object
    
    """
    
    with rasterio.Env():
        with rasterio.open(raster_path) as src:
            arr = src.read(masked=True)
            arr_ma, arr_ma_transform = mask(src, mask_feature, nodata=np.nan, crop=True,
                                              all_touched=False, invert=False)
    return arr_ma, arr_ma_transform

def ecef_normalizer(arr):
    """
    apply ecef normalization to 2d nd array
    
    For details can be found 
    https://www.statsmodels.org/stable/generated/statsmodels.distributions.empirical_distribution.ECDF.html
    
    ** need to import the model from statsmodel package
    from statsmodels.distributions.empirical_distribution import ECDF
    
    ** input to ECDF must be 1-dimensional
    
    Arguments
    arr: numpy nd array
        2-d nd array

    Returns
    arr_norm: normalized 2d nd array
    
    """
    
    # reshape to 1d
    arr_new = np.reshape(arr, -1)
    
    ecdf = ECDF(arr_new)
    arr_norm = ecdf(arr_new)
    
    # convert back to 2d nd array
    arr_norm = arr_norm.reshape(arr.shape[0], arr.shape[1])
    
    return arr_norm

def gaussian_filter_kernel5(arr):
    """
    apply gaussian filter to 2d nd array with a kernel size of 5x5
    
    **Need to import "gaussian_filter" module from scipy ndimage
    from scipy.ndimage import gaussian_filter
    
    For details can be found 
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter.html
    
    Parameters
    arr: numpy nd array
        2-d nd array

    Returns
    g_filtered: gaussian filtered (smoothed) 2d nd array
    
    """
    
    s = 4  # sigma (default=4)
    w = 5  # kernel size
    t = (((w - 1)/2)-0.5)/s # truncate

    g_filterd = gaussian_filter(arr, sigma=s, truncate=t)
    
    return g_filterd

def sieve_small_cluster(raster_path, out_path):

    """
    remove (sieve) small clusters using raterio (GDAL) sieve function

    ** threshold size for removing small area for production system is 8.573278464267728E-8 sqaure degrees
        that corresponds to 1072 square meters (~ 11 pixels at 10 m)
        but can be user defined with a value of number of pixels
    Argument
    raster_path: string
        raster file path

    Returns
    None

    """

    # Register GDAL and OGR drivers.
    with rasterio.Env():

        # Read a raster to be sieved.
        with rasterio.open(raster_path) as src:
            shade = src.read(1)
        mask = shade != np.nan
        
        # Sieve out features with 11 pixels or smaller.
        sieved = sieve(shade, 11, out=np.zeros(src.shape, src.dtypes[0]))

        # Print the number of shapes in the sieved raster.
        print("Sieved (11) shapes: %d" % len(list(shapes(sieved))))

        # Write out the sieved raster.
        kwargs = src.meta
        kwargs['transform'] = rasterio.transform.guard_transform(kwargs['transform'])

        with rasterio.open(out_path, 'w', **kwargs) as dst:
            dst.write(sieved, indexes=1)

def reproject_raster_latlon_to_utm_gdalwarp(infile_path, outfile_path, t_epsg, resolution):
    """
    reproject raster from lat/lon to UTM with a designated spatial resolution
    
    Arguments
    
    infile_path: string
        path to a input raster file
    outfile_path: string
        path to a input raster file
    t_epsg: string
        EPSG code in string format, for example epsg 32615 is UTM projection zone 15 in Norther Hemisphere
    resolution: string
        output resolution  in meter
    
    """
    cmd = "gdalwarp " + \
          "-of GTiff " + \
          "-tr " + resolution + ' ' + resolution + ' ' + \
          "-s_srs EPSG:4326 " + \
          "-t_srs EPSG:" + t_epsg + ' ' + \
          infile_path + ' ' + \
          outfile_path
    subprocess.check_call(cmd, shell=True)
    
def resample_raster_gdalwarp(infile_path, outfile_path, resolution):
    """
    Resample raster file. this is based on UTM projection which is the meter unit.
    The resolution must be the same unit to input raster. 
    
    Arguments
    infile_path: string
        path to input raster
    
    outfile_path: string
        path to output file
    
    resolution: string
        
    """
    cmd = "gdalwarp " + \
          "-of GTiff " + \
          "-tr " + resolution + ' ' + resolution + ' ' + \
          "-r bilinear " + \
          infile_path + ' ' + \
          outfile_path
    subprocess.check_call(cmd, shell=True)