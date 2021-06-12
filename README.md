# Zone_Vector

### This is a short project that support TCC Product in generating zone maps from a single raster file (VI or DEM) and converting it to vector format with ESRI shp file.

#### Main process includes 

1. read single raster file
2. normalize pixel values using the ECEF method
    * https://www.statsmodels.org/stable/generated/statsmodels.distributions.empirical_distribution.ECDF.html
    * reason for array normailzation is to make full stretch of pixel values within geometry
    * In some cases, pixel values might be homogeneous and not much variations
    * other normalization methods could be used as well
3. run k-means clustering algorithm and create zones
    * details of K-measn clustering by scikit learn; https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html
4. covert zone maps (in raster format) to vector file
    * a gdal library, gdal_polygonize is used and detail usage can be found here¶
https://gdal.org/programs/gdal_polygonize.html
5. (optional) removing small clusters
    * removes raster polygons smaller than a provided threshold size (in pixels) and replaces them with the pixel value of the largest neighbour polygon
    * threshold is **11 pixels**, but can be adjusted based on use cases
    * https://rasterio.readthedocs.io/en/latest/api/rasterio.features.html?highlight=sieve#rasterio.features.sieve
6. clipping zone maps to field geometry  

Mainworkflow can be found in the notebooks folder. 
* kmean_cluster_dem_vectorize_product.ipynb     
* custom_geospatial_utils.py contains a list of custom functions.     


There are 2 sample images with different spatial resolutions (pixel size). Coordinate Reference System for this example is the WGS84 and can be used other CRS system and input geotiff can be reprojected.
* us_3dep_10m_sample.tif - 10 m spatial resolution      
* ca_hidef_1m_sample.tif - 1 m spatial resolution       

Field geometry (in wkt format) file that corresponds to locations of each DEM file  
* sample_data_geometry.csv   

#### python library requirements
* rasterio
* gdal
* shapely
* KMeans
* ECEF
* subprocess
* custom_geospatial_tools (custom functions)