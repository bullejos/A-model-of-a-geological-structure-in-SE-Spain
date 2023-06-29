import numpy as np
import pandas as pd

import rasterio
import pyproj
def create_grid(xmin, ymin, xmax, ymax, xstep, ystep):
    """
    Create bounding coordinates for a regular 2d grid with numpy.
    These coordinates can then be used with shapely, pygeos and geopandas.
    """
    
    # Generate two evenly spaced ranges
    x = np.linspace(xmin, xmax, int(xstep+1))
    y = np.linspace(ymin, ymax, int(ystep+1))
    
    # Get center points for each range
    x = np.mean((x[1:], x[:-1]), axis=0)
    y = np.mean((y[1:], y[:-1]), axis=0)
    
    # Cartesian product of both ranges 
    x, y = np.meshgrid(x, y)
    x, y = x.ravel(), y.ravel()

    # Calculate delta x and delta y
    xd = (xmax - xmin) / (xstep * 2)
    yd = (ymax - ymin) / (ystep * 2)
    
    # Repeat center points and add respective delta
    grid = np.column_stack((x,y))
    grid = np.tile(grid, 2) + np.array([[-xd, -yd, xd, yd]])
    
    # Generate index
    grid_index = np.arange(0, xstep * ystep)
    
    # Output coordinates: xmin, ymin, xmax, ymax
    return grid, grid_index

def intersect_points_grid(xy, xmin, ymin, xmax, ymax, xstep, ystep, include_xmax=False, include_ymax=True):
    """
    Intersect 2d Points with a regular 2d Grid using NumPy.
    Returns the index of the Grid rectangle that each Point belongs to.
    """
    
    # Calculate delta x and delta y
    xd = (xmax - xmin) / xstep
    yd = (ymax - ymin) / ystep
    
    # Center points at zero
    xy = xy - np.array([[xmin, ymin]])
    
    # Get index of x and y coordinates
    xi = np.floor(xy[:,0] / xd) 
    yi = np.floor(xy[:,1] / yd)
    
    # Cast x and y index to 64 bit integer
    xi = xi.astype(np.int64)
    yi = yi.astype(np.int64)
    
    # Opened interval for x
    if include_xmax == True:
        xi = xi - (xy[:,0] % xd == 0)
    
    # Opened interval for y
    if include_ymax == True:
        yi = yi - (xy[:,1] % yd == 0)
        
    # Mask out-of-bound indexes on the x and y axis
    xi = np.ma.array(xi, mask=((xi < 0) | (xi >= xstep)))
    yi = np.ma.array(yi, mask=((yi < 0) | (yi >= ystep)))
        
    # Calculate index
    index = xi + (yi * xstep)
    
    # Unmask NumPy array
    index = index.filled(np.iinfo(np.int64).min)
    
    return index

class CopernicusDEM:
    """
    Extract elevation values for a Pandas dataframe of WGS84 coordinates,
    through raster files in Geotiff format provided by Copernicus programme.
    
    Since the rasters are of considerable size, 1000x1000km at 25 meter resolution,
    Rasterio's reading windows are used to process chunks individually, 
    limiting memory consuption.
    
    The process is completely implemented in Pandas, Numpy, Pyproj and Rasterio.
    
    - Source: Copernicus programme
    - Resource: European Digital Elevation Model (EU-DEM), version 1.1
    - URL: https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1
    """
    
    def __init__(self, raster_paths, grid_xstep=50, grid_ystep=50, elevation_alias='elevation'):
        """
        raster_paths: list of paths to the EU-DEM Geotiff raster files.
                      ['eu_dem_v11_E00N00.TIF', 'eu_dem_v11_E00N00.TIF'].
        
        grid_xstep: number of chunks in the x-axis used to create the reading windows. 
                    Must be diviseble by x-range. Affects performance and memory.
                    Default value is 50.
        
        grid_ystep: number of chunks in the y-axis used to create the reading windows. 
                    Must be diviseble by y-range. Affects performance and memory.
                    Default value is 50.
                    
        elevation_alias: name of output column with elevation values appended to the dataframe.         
        """
        
        # ETRS-LAEA grid parameters
        self.laea_xmin = 0
        self.laea_ymin = 0
        self.laea_xmax = 8000000
        self.laea_ymax = 6000000
        self.laea_xstep = 8
        self.laea_ystep = 6
        self.laea_include_xmax = False
        self.laea_include_ymax = True
        
        # User parameters and properties
        self.raster_paths = raster_paths
        self.grid_xstep = grid_xstep
        self.grid_ystep = grid_ystep
        
        # Input validations
        self._validate_grid_xstep()
        self._validate_grid_ystep()
        
        # Properties
        self.raster_crs = self._get_raster_crs()
        
        # Columns
        self.col_path = '_path'
        self.col_laea_tile = '_etrs_laea_tile'
        self.col_grid_index = '_grid_index'
        self.col_temp_key = '_temp_key'
        self.col_grid_xmin = '_grid_xmin'
        self.col_grid_ymin = '_grid_ymin'
        self.col_grid_xmax = '_grid_xmax'
        self.col_grid_ymax = '_grid_ymax'
        self.col_proj_x = '_x'
        self.col_proj_y = '_y'
        self.col_elevation = elevation_alias
                    
    def _validate_grid_xstep(self):
        """
        Check if xstep parameter is valid.
        """ 
        
        x_range = self.laea_xmax - self.laea_xmin
        x_step = self.laea_xstep * self.grid_xstep
        
        is_valid = x_range % x_step == 0
        assert_message = 'Parameter "grid_xstep" must be divisible by xrange ({}).'.format(x_range)
        
        assert is_valid, assert_message
        
    def _validate_grid_ystep(self):
        """
        Check if ystep parameter is valid.
        """ 
        
        y_range = self.laea_ymax - self.laea_ymin
        y_step = self.laea_ystep * self.grid_ystep
        
        is_valid = y_range % y_step == 0
        assert_message = 'Parameter "grid_ystep" must be divisible by yrange ({}).'.format(y_range)
        
        assert is_valid, assert_message
    
    def _get_raster_crs(self):
        """
        Read first geotiff raster file and retrieve Coordinate Reference System (CRS).
        """
        
        first_raster = self.raster_paths[0]
        
        with rasterio.open(first_raster) as src:
            crs = src.crs
        
        return crs
    
    def _get_etrs_laea_tile(self, df, xmin, ymin):
        """
        Convert rectangle bounds into ETRS-LAEA tile denomination ("E30N30").
        The first digit, at the kilometer level, of the xmin bound is East and of the ymin bound is North.
        """
        
        etrs_e = np.floor(df[xmin] / 10**6)
        etrs_e = (etrs_e * 10).astype(np.int64)
        etrs_e = 'E' + etrs_e.astype('str').str.zfill(2)
        
        etrs_n = np.floor(df[ymin] / 10**6)
        etrs_n = (etrs_n * 10).astype(np.int64)
        etrs_n = 'N' + etrs_n.astype('str').str.zfill(2)
        
        df[self.col_laea_tile] = etrs_e + etrs_n
        
        return df
        
    def _read_raster_metadata(self):
        """
        Read the metadata of every raster in the "raster_paths" parameter.
        Only the bounds of the raster and filename is kept.
        The bounds are then converted to the ETRS-LAEA tile denomination and dropped.
        """
        
        df = []
        for path in self.raster_paths:
            with rasterio.open(path) as src:
                bounds = dict(zip(('_raster_xmin', '_raster_ymin', '_raster_xmax', '_raster_ymax'), src.bounds))
                bounds[self.col_path] = path
                df.append(bounds)

        df = pd.DataFrame(df)
        df = self._get_etrs_laea_tile(df, xmin='_raster_xmin', ymin='_raster_ymin')
        df = df.drop(['_raster_xmin', '_raster_ymin', '_raster_xmax', '_raster_ymax'], axis=1)
        
        return df
    
    def _create_reading_grid(self):
        """
        Create the reading grid that will be the base for retrieving raster tiles.
        This grid must coincide and contain the original 1000x1000km ETRS-LAEA grid, and can 
        be subdivided into finer tiles, according to the "grid_xstep" and "grid_ystep" parameters.
        The ETRS-LAEA tile denomination is also added.
        Each tile of the created grid will have an index, represented by an intenger.
        """
        
        raster_grid, raster_grid_index = create_grid(xmin=self.laea_xmin, 
                                                     ymin=self.laea_ymin, 
                                                     xmax=self.laea_xmax, 
                                                     ymax=self.laea_ymax, 
                                                     xstep=self.laea_xstep*self.grid_xstep, 
                                                     ystep=self.laea_ystep*self.grid_ystep)
        
        df = pd.DataFrame(raster_grid, columns=[self.col_grid_xmin,self.col_grid_ymin,self.col_grid_xmax,self.col_grid_ymax])
        df[self.col_grid_index] = raster_grid_index
        df = self._get_etrs_laea_tile(df, xmin=self.col_grid_xmin, ymin=self.col_grid_ymin)
        
        return df
    
    def _get_grid_scope(self):
        """
        Grid scope is the intersection of the rasters included in input parameters and the reading grid.
        The ETRS-LAEA tile denomination is used as key in a database join.
        """
        
        df = self._create_reading_grid().merge(self._read_raster_metadata(), 'inner', self.col_laea_tile)
        
        df = df.drop(self.col_laea_tile, axis=1)
        
        return df
    
    def _project_points(self, df, lat_col, lon_col, input_crs):
        """
        Projects WGS84 coordinates into the raster CRS (EPSG:3035).
        The original coordinates are droped.
        """
        
        projection = pyproj.Transformer.from_crs(input_crs, self.raster_crs)
        
        df[[self.col_proj_x, self.col_proj_y]] = np.column_stack(projection.transform(df[lat_col].values, df[lon_col].values))
        
        df = df.drop([lat_col, lon_col], axis=1)
        
        return df
    
    def _get_grid_index(self, df):
        """
        Converts the projected points (x, y), into the respective reading grid index they belong to.
        """
        
        df[self.col_grid_index] = intersect_points_grid(xy=df[[self.col_proj_x, self.col_proj_y]].values, 
                                                        xmin=self.laea_xmin, 
                                                        ymin=self.laea_ymin, 
                                                        xmax=self.laea_xmax, 
                                                        ymax=self.laea_ymax, 
                                                        xstep=self.laea_xstep*self.grid_xstep, 
                                                        ystep=self.laea_ystep*self.grid_ystep,
                                                        include_xmax=self.laea_include_xmax,
                                                        include_ymax=self.laea_include_ymax)
        
        return df
    
    def _generate_temp_key(self, df):
        """
        Generates a temporary key that serves as reference 
        for the input dataset through the whole process.
        """
        
        df[self.col_temp_key] = range(len(df))
        
        return df
        
    def _get_elevation_pandas(self, df):
        """
        Function to be applied to partitions of the Pandas dataframe.
        It uses rasterio and all the predefined parameters to retrieve the elevation 
        of all points contained in the partition.
        The partitions are defined by grid index. This means that each window is only read once.
        The auxiliary columns that are no longer required are dropped to save memory.
        """
        
        path = df[self.col_path].iloc[0]
        xmin = df[self.col_grid_xmin].iloc[0]
        ymin = df[self.col_grid_ymin].iloc[0]
        xmax = df[self.col_grid_xmax].iloc[0]
        ymax = df[self.col_grid_ymax].iloc[0]
        
        cols_drop = [self.col_grid_index, 
                     self.col_path, 
                     self.col_grid_xmin, 
                     self.col_grid_ymin, 
                     self.col_grid_xmax, 
                     self.col_grid_ymax]
        
        df = df.drop(cols_drop, axis=1)
        
        with rasterio.open(path) as src:

            raster_affine = src.transform
            rasterio_window = rasterio.windows.from_bounds(left=xmin, 
                                                           bottom=ymin, 
                                                           right=xmax, 
                                                           top=ymax, 
                                                           transform=raster_affine)

            window_affine = src.window_transform(rasterio_window)
            index = rasterio.transform.rowcol(window_affine, df[self.col_proj_x], df[self.col_proj_y])
            df = df.drop([self.col_proj_x, self.col_proj_y], axis=1)
            
            elevation = src.read(1, window=rasterio_window)
            elevation = elevation[index]
            df[self.col_elevation] = elevation

            return df
    
    def _get_df_elevation(self, df, lat_col, lon_col, input_crs):
        """
        Base process to retrieve elevation for valid points.
        The reading grid with only tiles in scope is retrieved.
        The original df is copied, maintaining only the required columns (temp_id, lat, lon).
        Null latitude and longitude coordinates are filtered.
        The coordinates are then projected, and the grid index is obtained.
        The dataframe is filtered only to include points that fall within the bounds of the input rasters.
        The reading metadata, path and window bounds, are appended to the dataframe.
        Finaly the dataset is grouped or partioned through grid index and the elvation is obtained.
        """
        
        # Get grid tiles that are within scope of available input raster files.
        grid = self._get_grid_scope()
        
        # Copy original dataframe and retrieve necessary columns.
        # Drop null coordinates.
        df = df[[self.col_temp_key, lat_col, lon_col]].copy()
        df = df.dropna(subset=[lat_col, lon_col])
        
        # Project coordinates.
        # Get grid index.
        df = self._project_points(df, lat_col, lon_col, input_crs)
        df = self._get_grid_index(df)
        
        # Filter points that fall ouside of available grid bounds.
        # Join reading metadata.
        df = df.merge(grid[[self.col_grid_index]], 'inner', self.col_grid_index)
        df = df.merge(grid, 'inner', self.col_grid_index)
        
        # Group by grid index to obtain partitions.
        # Apply rasterio based function to get elevation.
        df = df.groupby(self.col_grid_index)
        df = df.apply(self._get_elevation_pandas)
        
        return df
    
    def get_elevation(self, df, lat_col, lon_col, input_crs='epsg:4326'):
        """
        User interface. Allows to retrieve elevation for a dataframe that contains WGS84 coordinates.
        The elevation column is appended to the input dataframe.
        All original columns of the dataframe are also returned (e.g. keys).
        Null input values values are supported and kept. 
        Elevation values for points that are out-of-bounds are treated as null.
        The raster uses a mask, a large negative value, around water bodies and the ocean.
        
        df: Pandas dataframe with coordinates and miscellaneous columns.
        lat_col: name of column that contains latitude.
        lon_col: name of columns that contains longitude.
        input_crs: input CRS, by default is WGS84 (EPSG:4326).
        """
        
        df = self._generate_temp_key(df)
        df = df.merge(self._get_df_elevation(df, lat_col, lon_col, input_crs), 'left', self.col_temp_key)
        df = df.drop(self.col_temp_key, axis=1)
        
        return df

