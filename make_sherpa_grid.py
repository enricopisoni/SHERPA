'''
Created on Jun 17, 2015

Create a grid of polygons in the database
The grid serves to calculate the surface of each cell and convert emission density into
absolute emissions

CREATE TABLE sherpa.grids
(
  id serial NOT NULL,
  grid_id integer,
  cell_id integer,
  the_cell geometry,
  CONSTRAINT pkey_id PRIMARY KEY (id)
)


select grid_id, count(cell_id) from sherpa.grids
group by grid_id
order by grid_id

query for surfaces

SELECT st_x(st_centroid(the_cell)) AS lon, st_y(st_centroid(the_cell)) AS lat, st_astext(the_cell), st_area(geography(the_cell))/1e6 as km2
FROM sherpa.grids
ORDER BY lon, lat

@author: degraba
'''

import psycopg2

host_name = "localhost"
database_name = "sherpa"
user_name = "postgres"
pw = "root"
schema_name = 'sherpa'
table_name = 'grids'


# database connection
try:
    conn = psycopg2.connect(host=host_name, database=database_name, user=user_name, password=pw);
except:
    print("I am unable to connect to the database")
cur = conn.cursor()

# truncate the tables
if 1 == 1:
    cur.execute('TRUNCATE TABLE sherpa.grids;')
    conn.commit()    


# grid properties
grid_id_lst = [1]
dlon = 0.125                # delta x in meters
dlat = 0.0625
lonmin = -10.5
lonmax = 37.5
latmin = 34.0    
latmax = 62.0    
srid_wgs84 = 4326

for grid_id in grid_id_lst:
    print('creating grid with dlon=%f and dlat=%f.' % (dlon, dlat))

    # create all cells of the grid from south to north and west to east
    # initialize cell id and coordinates
    cell_id = 1
    lon_west = lonmin
    lon_east = lonmin + dlon
    lat_south = latmin
    lat_north = latmin + dlat
    
    while lon_east <= lonmax:
        
        while lat_north <= latmax:
            
            # wkt for the polygon of the grid cell
            polygon_wkt = "ST_SetSRID(ST_GeomFromText('POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))'), %d)" \
            % (lon_west, lat_south, lon_east, lat_south, lon_east, lat_north, lon_west, lat_north, lon_west, lat_south, srid_wgs84)
            query = "INSERT INTO %s.%s(grid_id, cell_id, the_cell) VALUES (%d, %d, %s);" % (schema_name, table_name, grid_id, cell_id, polygon_wkt)
            # print(query)
            cur.execute(query)
            conn.commit()
            cell_id += 1
            # go one cell north
            lat_south += dlat
            lat_north += dlat
    
        # go one cell east an go back south
        lon_east += dlon
        lon_west += dlon
        lat_south = latmin
        lat_north = latmin + dlat
    
    # create table ifdm_ia.emission_gridding as select ST_Intersection(the_road, the_cell) from ifdm.lijnbroninvoer_antwerp, (select * from ifdm_ia.grids where cell_id=119) g