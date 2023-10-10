"""
Fast watershed delineation for anywhere on the Earth's surface.
Creates watersheds for outlet points in an input CSV file.

Interactive demo: https://mghydro.com/watersheds

Usage: Carefully edit the file config.py, then run delineator.py

Capable of outputting to many formats, and includes a handy viewer for web browsers.

In "low-resolution" mode, creates watersheds that include some area downstream of the outlet, because
of the way the watershed is created by merging pre-existing unit catchments, from
MERIT-Hydro vector watersheds GIS data from Princeton Univ. The errors average
about 20 km², but can occasionally be much larger. For large watersheds, these errors may not be
important.

Use "high-resolution" mode for greater accuracy, especially for smaller watersheds. This requires
more input data and takes more time. To use detailed mode, you must first follow the directions
to download MERIT-Hydro raster data from U. of Tokyo.

For comments or questions, please contact the author: Matthew Heberger, matt@mghydro.com
or create an Issue on the GitHub repo: https://github.com/mheberger/delineator
"""
import warnings

import numpy as np
import pickle

warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import geopandas as gpd
import os
import re
from shapely.geometry import Point, Polygon, box
import shapely.ops
from shapely.wkt import loads
import sigfig  # for formatting numbers to significant digits
from py.fast_dissolve import dissolve_geopandas, fill_geopandas
import pyproj
from functools import partial
from config import *
from py.mapper import make_map, create_folder_if_not_exists

if PLOTS:
    import matplotlib.pyplot as plt

if HIGH_RES:
    import py.merit_detailed


gpd.options.use_pygeos = True

# The WGS84 projection string, used in a few places
PROJ_WGS84 = 'EPSG:4326'


def validate(gages_df: pd.DataFrame) -> bool:
    """
    After we have read in the gages CSV, check whether it appears to be OK
    if so, return True

    """
    cols = gages_df.columns
    required_cols = ['id', 'lat', 'lng']
    for col in required_cols:
        if col not in cols:
            raise Exception(f"Missing column in CSV file: {col}")


    # Check that the ids are all unique
    if len(gages_df['id'].unique()) != len(gages_df):
        raise Exception("Each id in your CSV file must be unique.")

    # Check that lat, lng are numeric
    fields = ['lat', 'lng']
    for field in fields:
        if gages_df[field].dtype != 'float64':
            raise Exception(f"In outlets CSV, the column {field} is not numeric.")

    # Check that all the lats are in the right range
    lats = gages_df["lat"].tolist()
    lngs = gages_df["lng"].tolist()

    if not all(lat > -60 for lat in lats):
        raise Exception("All latitudes must be greater than -60°")

    if not all(lat < 85 for lat in lats):
        raise Exception("All latitudes must be less than 85°")

    if not all(lng > -180 for lng in lngs):
        raise Exception("All longitudes must be greater than -180°")

    if not all(lng < 180 for lng in lngs):
        raise Exception("All longitudes must be less than 180°")

    # Check that every row has an id
    ids = gages_df["id"].tolist()

    if not all(len(str(wid)) > 0 for wid in ids):
        raise Exception("Every watershed outlet must have an id in the CSV file")

    return True


def validate_search_distance():
    """
    For the pour point relocation routine. SEARCH_DIST governs how far away we'll look for
    a unit catchment with an upstream area close to the a priori estimate.
    This function just checks whether this is valid.
    Does not return anything, just throws an error if it is out of range.
    """
    if not isinstance(SEARCH_DIST, int) and not isinstance(SEARCH_DIST, float):
       raise Exception(f"SEARCH_DIST must be a number. We got {SEARCH_DIST}")

    if SEARCH_DIST < 0.0:
        raise Exception("SEARCH distance in config.py must be a positive number.")

    if SEARCH_DIST > 0.25:
        raise Exception("SEARCH_DIST is unrealistically high. It should be in decimal degrees, and must be less "
                        "than 0.25. In config.py, you entered {SEARCH_DIST}")


def get_area(poly: Polygon) -> float:
    """
    Projects a Shapely polygon in raw lat, lng coordinates, and calculates its area
    Argss:
        poly: Shapely polygon
    :return: area in km²
    """
    projected_poly = shapely.ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(init='EPSG:4326'),
            pyproj.Proj(
                proj='aea',
                lat_1=poly.bounds[1],
                lat_2=poly.bounds[3]
            )
        ),
        poly)

    # Get the area in m^2
    return projected_poly.area / 1e6


def delineate():
    """
    THIS is the Main watershed delineation routine
    Make sure to set the variables in `config.py` before running.

    Reads a list of outlet points from a .csv file,
    then finds their watersheds or drainage basins,
    using hybrid of vector- and raster-based methods.

    Outputs geodata (.shp, .gpkg, etc.) and optionally a CSV file with a summary of results
    (shows the watershed id, names, and areas).

    Optionally creates an HTML page with a handy map viewer to review the results.
    """

    # Regular expression used to find numbers so I can round lat, lng coordinates in GeoJSON files to make them smaller
    simpledec = re.compile(r"\d*\.\d+")

    def mround(match):
        # For rounding the coordinates in GeoJSON files to make them smaller
        return "{:.5f}".format(float(match.group()))

    def addnode(B, node):
        """"
        Recursive function to assemble the list of upstream unit catchments
        B is a list of the unit catchments that make up our watershed.
        List items are `comid`s, unique identifiers of each unit catchment.
        """
        # first, append the node to the basin
        B.append(node)

        # next, check whether the fields up1, up2, up3, and up4 contain a node ID
        up1 = rivers_gdf['up1'].loc[node]
        if up1 != 0:
            addnode(B, up1)

        up2 = rivers_gdf['up2'].loc[node]
        if up2 != 0:
            addnode(B, up2)

        up3 = rivers_gdf['up3'].loc[node]
        if up3 != 0:
            addnode(B, up3)

        up4 = rivers_gdf['up4'].loc[node]
        if up4 != 0:
            addnode(B, up4)

    def find_close_catchment():
        """
        Part of my simple pour point relocation method. If the outlet falls in a unit catchment
        whose area is a mismatch of our a priori estimate of the upstream area,
        look around the neighborhood for a unit catchment whose area matches more closely
        """
        # Find the river segment that is within a bounding box around the point
        dist = 0.01
        max_dist = MAX_DIST  # How far away can we look before we give up?

        # Keep track of how many river segments were found with each increase in distance.
        num_segments_found = [0]

        iteration = 1
        while True:
            find_box = box(lng - dist, lat - dist, lng + dist, lat + dist)
            possible_matches_index = list(rivers_gdf.sindex.intersection(find_box.bounds))
            possible_matches = rivers_gdf.iloc[possible_matches_index]
            precise_matches = possible_matches[possible_matches.intersects(find_box)]
            segments_found = precise_matches
            num_segments_found.append(len(segments_found.index))

            # If any new segments were found on this iteration...
            if num_segments_found[iteration] > num_segments_found[iteration - 1]:
                # Calculate the percent difference in the area for the found river segments and the gage
                segments_found['pd'] = ((segments_found['uparea'] - area_reported) / area_reported).abs()

                # Get some info on the "best matching" river segment
                min_pd = segments_found['pd'].min()

                # Check whether they are "good enough"
                if min_pd < AREA_MATCHING_THRESHOLD:
                    COMID = segments_found['pd'].idxmin()
                    uparea = round(segments_found['uparea'][COMID], 0)
                    if VERBOSE:
                        print("  (X) Found a river reach at a distance of {}, "
                              "area difference: {:,.0f}%".format(dist, min_pd * 100))
                    return COMID, uparea

            if dist > max_dist:
                # If we've gone out a certain radius around the gage, and still haven't found a river
                # segment with a closely matching upstream area, raise some kind of error message.
                print('  (!) Could not find a river segment with closely matching upstr.'
                      ' area within %s degrees of gage' % round(dist, 2))

                return None, None
            else:
                # Otherwise, expand the search radius and keep looking
                dist += 0.01
                iteration += 1

    def plot_basins(suffix: str):
        """
        Makes a plot of the unit catchments that are in the watershed

        It is all the upstream unit catchments, and the *split* terminal unit catchment.

        """
        # subbasins_gdf.plot(column='area', edgecolor='gray', legend=True)
        [fig, ax] = plt.subplots(1, 1, figsize=(10, 8))

        # Plot each unit catchment with a different color
        for x in subbasins_gdf.index:
            color = np.random.rand(3, )
            subbasins_gdf.loc[[x]].plot(facecolor=color, edgecolor=color, alpha=0.5, ax=ax)

        # Plot the gage point
        gages_joined.iloc[[i]].plot(ax=ax, c='red', edgecolors='black')

        if suffix == "post" and HIGH_RES:
            plt.scatter(x=lng_snap, y=lat_snap, c='cyan', edgecolors='black')
            plt.title(f"Showing the {len(subbasins_gdf)-1} upstream unit catchments and split terminal unit catchment")
        else:
            plt.title(f"Found {len(subbasins_gdf)} unit catchments for watershed id = {wid}")

        plt.savefig(f"plots/{wid}_vector_unit_catchments_{suffix}.png")
        plt.close(fig)

    # Check that the OUTPUT directories are there. If not, try to create them.
    folder_exists =  create_folder_if_not_exists(OUTPUT_DIR)
    if not folder_exists:
        raise Exception(f"No folder for output. Stopping")

    # Check for the folder to put Python PICKLE files
    if PICKLE_DIR != "":
        folder_exists =  create_folder_if_not_exists(OUTPUT_DIR)
        if not folder_exists:
            raise Exception(f"No folder for pickle files. Stopping")

    # Check if the MAP_FOLDER is there
    if MAKE_MAP:
        folder_exists = create_folder_if_not_exists(MAP_FOLDER)
        if not folder_exists:
            raise Exception(f"No folder for the map files. Stopping")

    # Check that the CSV file is there
    if not os.path.isfile(OUTLETS_CSV):
        raise Exception(f"Could not your outlets file at: {OUTLETS_CSV}")

    # Read the outlet points CSV file and put into a Pandas DataFrame
    # (I call the outlet points gages, because I usually in delineated watersheds at streamflow gages)
    if VERBOSE: print(f"Reading your outlets data in: {OUTLETS_CSV}")
    gages_df = pd.read_csv(OUTLETS_CSV, header=0, dtype={'id': 'str', 'lat': 'float', 'lng': 'float'})

    # Check that the CSV file includes at a minimum: id, lat, lng and that all values are appropriate
    validate(gages_df)

    # Get the number of points, for status messages
    n_gages = len(gages_df)

    # Boolean vars to track whether user's outlets file had fields `area` and/or `name`
    bAreas = 'area' in gages_df
    bNames = 'name' in gages_df

    # Add some extra fields to the gages dataframe
    if bAreas:
        gages_df['area_reported'] = pd.to_numeric(gages_df['area'])
        gages_df.drop(['area'], axis=1, inplace=True)

    # If we are doing detailed delineation with raster data, we'll keep track of the "snapped" pour point
    # pysheds will always move the point a little bit, so that it coincides with the gridded data.
    if HIGH_RES:
        gages_df['lat_snap'] = np.nan
        gages_df['lng_snap'] = np.nan
        gages_df['snap_dist'] = 0

    # Convert gages_df to a GeoPandas GeoDataFrame (adds geography, lets us do geo. operations)
    coordinates = [Point(xy) for xy in zip(gages_df['lng'], gages_df['lat'])]
    points_gdf = gpd.GeoDataFrame(gages_df, crs=PROJ_WGS84, geometry=coordinates)
    # The line above has the surpsising side effect of adding a geometry column to gages_df (!)
    # Since we don't want or need this, drop the column.
    # No longer needed with GeoPandas v 0.14
    #gages_df.drop(['geometry'], axis=1, inplace=True)

    if VERBOSE: print("Finding out which Level 2 megabasin(s) your points are in")
    # This file has the merged "megabasins" in it
    merit_basins_shp = 'data/shp/basins_level2/merit_hydro_vect_level2.shp'
    megabasins = gpd.read_file(merit_basins_shp)
    # The CRS string in the shapefile is EPSG 4326 but does not match verbatim
    megabasins.to_crs(PROJ_WGS84, inplace=True)
    if not megabasins.loc[0].BASIN == 11:
        raise Exception("An error occurred loading the Level 2 basins shapefile")

    # Overlay the gage points on the Level 2 Basins polygons to find out which
    # PFAF_2 basin each point falls inside of, using a spatial join
    if SEARCH_DIST == 0:
        gages_basins_join = gpd.sjoin(points_gdf, megabasins, op="within")
    else:
        # Better results obtained with a "nearest shape"
        # This line generates a warning about how its bad to use distances in unprojected geodata. OK
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=UserWarning)
            gages_basins_join = gpd.sjoin_nearest(points_gdf, megabasins, max_distance=SEARCH_DIST)

    # Needed to set this option in order to avoid a warning message in Geopandas.
    # https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
    pd.options.mode.chained_assignment = None  # default='warn'

    # Keep a record of gages for which we could not find a matching river segment.
    # Dict failed: key = id, value = string, explanation of failure
    failed = {}

    # Get a list of the DISTINCT Level 2 basins, and a count of how many gages in each.
    basins_df = gages_basins_join.groupby("BASIN").id.nunique()
    basins = basins_df.index.tolist()

    if VERBOSE: print(f"Your watershed outlets are in {len(basins)} basin(s)")

    # Find any outlet points that are not in any Level 2 basin, and add these to the fail list
    # Look for any rows that are in gages_df that are not in basins_df
    matched_ids = gages_basins_join['id'].tolist()
    ids = gages_df['id'].tolist()
    for wid in ids:
        if wid not in matched_ids:
            failed[wid] = "Point not located in any Level 2 basin."

    # Now add fields to gages_df so we can reuse it to create a table to output to CSV
    gages_df.set_index('id', inplace=True)
    gages_df['area_calc'] = 0
    gages_df['result'] = "failed"
    if bAreas:
        gages_df['perc_diff'] = 0

    gages_counter = 0

    # Iterate over the basins so that we only have
    # to open up each Level 2 Basin shapefile once, and handle all of the gages in it
    for basin in basins:
        # Create a dataframe of the gages_basins_join in that basins
        gages_in_basin = gages_basins_join[gages_basins_join["BASIN"] == basin]
        num_gages_in_basin = len(gages_in_basin)
        if VERBOSE: print("\nBeginning delineation for %s outlet point(s) in Level 2 Basin #%s." % (num_gages_in_basin, basin))

        if HIGH_RES:
            catchments_gdf = load_gdf("catchments", basin, True)
            catchments_lowres_gdf = None
        else:
            catchments_gdf = load_gdf("catchments", basin, False)

        # The network data is in the RIVERS file rather than the CATCHMENTS file
        # (this is just how the MeritBASIS authors did it)
        if VERBOSE: print('Reading data table for rivers in basin %s' % basin)
        rivers_gdf = load_gdf("rivers", basin, True)

        # Perform a Spatial join on gages (points) and unit catchments (polygons)
        # to find the corresponding unit catchment for each gage
        # Adds the fields COMID and unitarea
        if VERBOSE: print(f"Performing spatial join on {num_gages_in_basin} outlet points in basin #{basin}")
        gages_in_basin.drop(['index_right'], axis=1, inplace=True)
        validate_search_distance()
        if SEARCH_DIST == 0:
            gages_joined = gpd.sjoin(gages_in_basin, catchments_gdf, how="inner", predicate="intersects")
        else:
            # This line generates a warning about how its bad to use distances in unprojected geodata. OK
            with warnings.catch_warnings():
                warnings.simplefilter(action='ignore', category=UserWarning)
                gages_joined = gpd.sjoin_nearest(gages_in_basin, catchments_gdf, max_distance=SEARCH_DIST)

            gages_joined.rename(columns={"index_right": "COMID"}, inplace=True)

        # For any gages for which we could not find a unit catchment, add them to failed
        gages_matched = gages_joined['id'].tolist()
        gage_basin_ids = gages_in_basin['id'].tolist()
        for wid in gage_basin_ids:
            if wid not in gages_matched:
                failed[wid] = f"Could not assign to a unit catchment in Level 2 basin #{basin}"

        # Revise the number of gages in the basin, based on those which have a matching COMID
        num_gages_in_basin = len(gages_joined)

        # Iterate over the gages and assemble the watershed
        for i in range(0, num_gages_in_basin):
            gages_counter += 1
            if VERBOSE: print(f"\n* Delineating watershed {gages_counter} of {n_gages}, with outlet id = {wid}")

            # Reset the local boolean flag for high-res mode. If the watershed is too big, script will switch to low.
            bool_high_res = HIGH_RES

            # Let wid be the watershed ID. Get the lat, lng coords of the gage.
            wid = gages_joined['id'].iloc[i]
            lat = gages_joined['lat'].iloc[i]
            lng = gages_joined['lng'].iloc[i]

            # The terminal comid is the unit catchment that contains (overlaps) the outlet point
            terminal_comid = gages_joined['COMID'].iloc[i]

            # Get the upstream area of the unit catchment we found, according to MERIT-Basins
            up_area = rivers_gdf.loc[terminal_comid].uparea

            # If MATCH_AREAS is True and the user provided the area in the outlets CSV file,
            # the script will check whether the upstream area of the unit catchment is a good match.
            # If the areas do not match well, look around the neighborhood for another unit catchment
            # whose area is a closer match to what we think it is.
            if bAreas and MATCH_AREAS:
                area_reported = gages_df.loc[wid].area_reported
                PD_area = abs((area_reported - up_area) / area_reported)
                if PD_area > AREA_MATCHING_THRESHOLD:
                    if VERBOSE:
                        print("Outlet point is in a unit catchment whose area is not a close match.")
                        print("Searching neighborhood for a river reach with a more closely matching upstream area")
                    candidate_comid, up_area = find_close_catchment()
                    if candidate_comid is None:
                        failed[wid] = "Could not find a nearby river reach whose upstream area is " \
                                      "within {}% of reported area of {:,.0f} km²" \
                            .format(AREA_MATCHING_THRESHOLD * 100, area_reported)
                        continue
                    else:
                        terminal_comid = candidate_comid

            # Let B be the list of unit catchments (and river reaches) that are in the basin
            B = []

            # Add the first node, and the rest will be added recursively
            addnode(B, terminal_comid)
            if VERBOSE: print(f"  found {len(B)} unit catchments in the watershed")

            # If the watershed is too big, revert to low-precision mode.
            if HIGH_RES and up_area > LOW_RES_THRESHOLD:
                if VERBOSE: print(f"Watershed for id = {wid} is larger than LOW_RES_THRESHOLD = {LOW_RES_THRESHOLD}. "
                                  "SWITCHING TO LOW-RESOLUTION MODE.")
                bool_high_res = False
                # If we just flipped to low-res mode, check if the low-res unit catchment polygons are loaded.
                if catchments_lowres_gdf is None:
                    catchments_lowres_gdf = load_gdf("catchments", basin, False)

                subbasins_gdf = catchments_lowres_gdf.loc[B]
            else:
                # Create a new geodataframe containing only the unit catchments_gdf that are in the list B
                # In high precision mode, we will update the geometry of the terminal unit catchment.
                subbasins_gdf = catchments_gdf.loc[B]

            # Make a plot of the selected unit catchments
            if PLOTS:
                plot_basins("pre")

            # In detailed mode,
            if bool_high_res:
                if VERBOSE: print("Performing detailed raster-based delineation for "
                                  "the downstream portion of the watershed")
                # Let split_catchment_poly be the polygon of the terminal unit catchment
                assert terminal_comid == B[0]
                catchment_poly = subbasins_gdf.loc[terminal_comid].geometry
                bSingleCatchment = len(B) == 1
                split_catchment_poly, lat_snap, lng_snap = py.merit_detailed.split_catchment(wid, basin, lat, lng,
                                                                                             catchment_poly,
                                                                                             bSingleCatchment)
                if split_catchment_poly is None:
                    failed[wid] = "An error occured in pysheds detailed delineation."
                    continue
                else:
                    # Create a temporary GeoDataFrame to create the geometry of the split catchment polygon
                    # This is just a shortcut method to transfer it to our watershed's subbasins GDF
                    split_gdf = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[split_catchment_poly])
                    split_geom = split_gdf.loc[0, 'geometry']
                    subbasins_gdf.loc[terminal_comid, 'geometry'] = split_geom

            if PLOTS:
                plot_basins("post")

            if VERBOSE: print("Dissolving...")
            # mybasin_gs is a GeoPandas GeoSeries
            mybasin_gs = dissolve_geopandas(subbasins_gdf)

            if FILL:
                # Fill donut holes in the watershed polygon
                # Recall we asked the user for the fill threshold in terms of number of pixels
                PIXEL_AREA = 0.000000695  # Constant for the area of a single pixel in MERIT-Hydro, in decimal degrees
                area_max = FILL_THRESHOLD * PIXEL_AREA
                mybasin_gs = fill_geopandas(mybasin_gs, area_max=area_max)

            if SIMPLIFY:
                # Simplify the geometry. GeoPandas uses the simple Douglas-Peuker algorithm
                mybasin_gs = mybasin_gs.simplify(tolerance=SIMPLIFY_TOLERANCE)

            # Let mybasin_gdf be a GeoPandas DataFrame with the geometry, and the id and area of our watershed
            mybasin_gdf = gpd.GeoDataFrame(geometry=mybasin_gs)
            mybasin_gdf['id'] = wid
            basin_poly = mybasin_gdf.geometry.values[0]
            up_area = get_area(basin_poly)
            # If the user gave a name and an a priori area to the watershed, include it in the output
            if bNames:
                mybasin_gdf['name'] = gages_df.loc[wid, 'name']
            if bool_high_res:
                mybasin_gdf['result'] = "High Res"
                gages_df.at[wid, 'result'] = "high res"

            else:
                mybasin_gdf['result'] = "Low Res"
                gages_df.at[wid, 'result'] = "low res"

                snapped_outlet = rivers_gdf.loc[terminal_comid].geometry.coords[0]
                lat_snap = snapped_outlet[1]
                lng_snap = snapped_outlet[0]

            # Get the (approx.) snap distance
            geod = pyproj.Geod(ellps='WGS84')
            snap_dist = geod.inv(lng, lat, lng_snap, lat_snap)[2]
            gages_df.at[wid, 'snap_dist'] = sigfig.round(snap_dist, 2)
            gages_df.at[wid, 'lat_snap'] = round(lat_snap, 3)
            gages_df.at[wid, 'lng_snap'] = round(lng_snap, 3)


            up_area = sigfig.round(up_area, 3)
            mybasin_gdf['area_calc'] = up_area
            gages_df.at[wid, 'area_calc'] = up_area

            if bAreas:
                area_reported = gages_df.loc[wid].area_reported
                mybasin_gdf['area_reported'] = area_reported
                perc_diff = sigfig.round((up_area - area_reported) / area_reported * 100, 2)
                gages_df.at[wid, 'perc_diff'] = perc_diff

            # SAVE the Watershed to disk as a GeoJSON file or a shapefile
            if VERBOSE: print(f' Writing output for watershed {wid}')
            outfile = f"{OUTPUT_DIR}/{wid}.{OUTPUT_EXT}"

            # This line rounds all the vertices to fewer digits. For text-like formats GeoJSON or KML, makes smaller
            # files with minimal loss of precision. For other formats (shp, gpkg), doesn't make a difference in file size
            if OUTPUT_EXT.lower() in ['geojson', 'kml']:
                mybasin_gdf.geometry = mybasin_gdf.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

            if OUTPUT_EXT != "":
                with warnings.catch_warnings():
                    warnings.simplefilter(action='ignore', category=UserWarning)
                    mybasin_gdf.to_file(outfile)

            # Create the HTML Viewer Map?
            # We have to write a second, slightly different version of the GeoJSON files,
            # because we need it in a .js file assigned to a variable, to avoid cross-origin restrictions
            # of modern web browsers.
            if MAKE_MAP:
                watershed_js = f"{MAP_FOLDER}/{wid}.js"
                with open(watershed_js, 'w') as f:
                    s = f"gage_coords = [{lat}, {lng}];\n"
                    f.write(s)
                    s = f"snapped_coords = [{lat_snap}, {lng_snap}];\n"
                    f.write(s)

                    f.write("basin = ")
                    f.write(mybasin_gdf.to_json())

                if MAP_RIVERS:
                    myrivers_gdf = rivers_gdf.loc[B]

                    # Keep only the fields lengthkm and order
                    myrivers_gdf = myrivers_gdf[['lengthkm', 'order', 'geometry']]

                    # Filter out the little headwater streams in large watersheds.
                    max_order = myrivers_gdf.order.max()
                    min_order = max_order - NUM_STREAM_ORDERS
                    # Drop rows where order < min_order
                    myrivers_gdf = myrivers_gdf[myrivers_gdf.order >= min_order]
                    myrivers_gdf = myrivers_gdf.round(1)
                    myrivers_gdf.geometry = myrivers_gdf.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))
                    rivers_js = f"{MAP_FOLDER}/{wid}_rivers.js"
                    with open(rivers_js, 'w') as f:
                        f.write("rivers = ")
                        f.write(myrivers_gdf.to_json())

    # CREATE OUTPUT.CSV, a data table of the outputs
    # id, status (hi, low, failed), name, area_reported, area_calculated
    if OUTPUT_CSV:
        output_csv_filename = f"{OUTPUT_DIR}/OUTPUT.csv"
        gages_df.to_csv(output_csv_filename)

    # FAILED.csv: If there were any failures, write this to a separate CSV file
    if len(failed) > 0:
        print(f"### FAILED to find watersheds for {len(failed)} locations. Check FAILED.csv for info.")

        failfile = f"{OUTPUT_DIR}/FAILED.csv"
        with open(failfile, 'w') as f:
            f.write("ID, EXPLANATION\n")
            for k, v in failed.items():
                f.write(f'{k},"{v}"\n')

    # If the user wants the browser map, make it
    if MAKE_MAP:
        if VERBOSE: print("* Creating viewer.html *")
        make_map(gages_df)

    # Finished, print a little status message
    if VERBOSE: print(f"It's over! See results in {output_csv_filename}")


def get_pickle_filename(geotype: str, basin: int, high_resolution: bool) -> str:
    """Simple function to get the standard filename for the pickle files used by this project.
    The filenames look like this:
       PICKLE_DIR/catchments_##_hires.pkl
       PICKLE_DIR/catchments_##_lores.pkl

       PICKLE_DIR/rivers_##_hires.pkl
       PICKLE_DIR/rivers_##_lores.pkl

    where ## is the megabasin number (11-91)

    """

    if high_resolution:
        resolution_str = 'hires'
    else:
        resolution_str = 'lores'
    fname = f'{PICKLE_DIR}/{geotype}_{basin}_{resolution_str}.pkl'
    return fname


def load_gdf(geotype: str, basin: int, high_resolution: bool) -> gpd.GeoDataFrame:
    """
    Returns the unit catchments vector polygon dataset as a GeoDataFrame
    Gets the data from the MERIT-Basins shapefile the first time,
    and after that from a saved .pkl file on disk.
    Uses some global parameters from config.py

    :param geotype: either "catchments" or "rivers" depending on which one we want to open.
    :param basin: the Pfafstetter level 2 megabasin, an integer from 11 to 91
    :param high_resolution: True to load the standard (high-resolution) file,
      False to load the low-resolution version (for faster processing, slightly less accurate results)

    :return: a GeoPandas GeoDataFrame

    """

    # First, check for the presence of a pickle file
    if PICKLE_DIR != '':
        pickle_fname = get_pickle_filename(geotype, basin, high_resolution)
        if os.path.isfile(pickle_fname):
            if VERBOSE: print(f"Fetching BASIN # {basin} catchment data from pickle file.")
            gdf = pickle.load(open(pickle_fname, "rb"))
            return gdf

    # Open the shapefile for the basin
    if geotype == "catchments":
        if high_resolution:
            directory = HIGHRES_CATCHMENTS_DIR
        else:
            directory = LOWRES_CATCHMENTS_DIR
        shapefile = f"{directory}/cat_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01.shp"
    elif geotype == "rivers":
        shapefile = f"{RIVERS_DIR}/riv_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01.shp"

    if not os.path.isfile(shapefile):
        raise Exception(f"Could not find the file: {shapefile}")

    if VERBOSE: print(f"Reading geodata in {shapefile}")
    gdf = gpd.read_file(shapefile)
    gdf.set_index('COMID', inplace=True)

    # This line is necessary because some of the shapefiles provided by reachhydro.com do not include .prj files
    gdf.set_crs(PROJ_WGS84, inplace=True, allow_override=True)

    # Before we exit, save the GeoDataFrame as a pickle file, for future speedups!
    save_pickle(geotype, gdf, basin, high_resolution)
    return gdf


def save_pickle(geotype: str, gdf: gpd.GeoDataFrame, basin: int, high_resolution: bool):
    # If we loaded the catchments from a shapefile, save the gdf to a pickle file for future speedup
    if PICKLE_DIR != '':

        # Check whether the GDF has a spatial index.
        # Note: I don't think this is ever necessary. Since version 0.7.0 (March 2020), GeoPandas
        # creates a spatial index by default.
        has_spatial_index = hasattr(gdf, 'sindex')
        if not has_spatial_index:
            gdf.sindex.create_index()

        # Get the standard project filename for the pickle files.
        pickle_fname = get_pickle_filename(geotype, basin, high_resolution)
        if not os.path.isfile(pickle_fname):
            if VERBOSE: print(f"Saving GeoDataFrame to pickle file: {pickle_fname}")
            try:
                pickle.dump(gdf, open(pickle_fname, "wb"))
            except:
                raise Warning("Could not save pickle file to: {pickle_fname}")


if __name__ == "__main__":
    delineate()
