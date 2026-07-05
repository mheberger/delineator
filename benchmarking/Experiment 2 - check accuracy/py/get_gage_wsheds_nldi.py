# Get the watershed boundaries for gages from the NLDI API and save them
# as shapefiles.

import pandas
import geopandas
import requests
import matplotlib.pyplot as plt

# First, load the list of gages from a CSV file.
# We'll put this into a Pandas dataframe for convenience
csv_file = r"C:\Users\mheberger\Documents\delineator\PAPER\usgs_gage_sampling\miss.csv"

#gages_df = pd.read_csv(csv_file)
#print(len(gages_df))

with open(csv_file, 'r') as f:
    f.readline()  # We don't need line #1
    lines = f.readlines()

gages = []
names = []
for line in lines:
    data = line.split(",")
    gage_id = data[0]
    gages.append(gage_id)
    name = data[1]
    names.append(name)

print(len(gages))

# Now, get the watershed GeoJSON from the NLDI API.
i = 0
for gage in gages:
    print("Processing {}".format(i))
    url = "https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-{}/basin".format(gage)
    r = requests.get(url=url)
    #data = r.json()
    watershed_gdf = geopandas.read_file(r.text)

    #print(watershed_gdf.head())

    #ax = watershed_gdf.plot(alpha=0.5)

    #plt.title("Watershed boundary")
    #plt.show()

    fname = "C:/Users/mheberger/Documents/delineator/PAPER/miss_vick.gpkg".format(gage)
    watershed_gdf.to_file(fname)

    i += 1


