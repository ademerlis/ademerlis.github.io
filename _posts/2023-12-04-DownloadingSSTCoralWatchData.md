---
layout: post
title: Downloading SST Data from CoralWatch
date: '2023-12-04'
categories: Coding
tags: [Coding]
---

Asking ChatGPT to define everything I need: 

"Degree Heating Weeks (DHW) is a measure used to assess the thermal stress on coral reefs, which is a significant factor in coral bleaching events. It's calculated by accumulating the amount of time sea surface temperatures exceed the long-term historical maximum mean temperature for a particular location. Here's how you calculate DHW:

1. **Determine the Historical Maximum Mean Temperature (MMM):** This is the highest average monthly sea surface temperature for a given location over a long-term period (usually over several decades). This temperature serves as a threshold.

2. **Monitor Daily Sea Surface Temperatures (SST):** Obtain daily SST data for your location. This data is often available from satellite observations or oceanographic institutions.

3. **Calculate the Temperature Anomaly:** For each day, determine how much the SST exceeds the MMM. This is the temperature anomaly. If the SST is lower than the MMM, the anomaly is considered zero for that day.

4. **Accumulate the Anomalies Over Time:** Degree Heating Weeks are calculated over a rolling 12-week (84-day) period. For each day, add the temperature anomaly to a cumulative sum. Continue this process for each day in the 12-week period.

5. **Express the Total as Weeks:** The total accumulated temperature anomaly is then divided by 7 to express it in terms of "degree heating weeks." This unit represents the number of weeks the sea temperature has been above the MMM.

6. **Update Regularly:** DHW is a rolling measure, so it's updated regularly (daily or weekly) to reflect ongoing temperature conditions.

DHW is a critical tool for marine scientists and environmental managers in monitoring and responding to the impacts of climate change on coral reef ecosystems. High DHW values indicate prolonged thermal stress, which can lead to coral bleaching and other ecological impacts.

Determining the Historical Maximum Mean Temperature (MMM) for a specific location involves analyzing long-term sea surface temperature data. The MMM is the highest average monthly sea surface temperature typically calculated over a baseline period of several decades. Here's a step-by-step process to determine the MMM:

1. **Collect Long-term Sea Surface Temperature Data:** Obtain historical sea surface temperature (SST) data for your location. This data is often available from meteorological organizations, oceanographic institutes, or satellite observations. The dataset should ideally span at least 20-30 years to capture climate variability.

2. **Select the Baseline Period:** Choose a baseline period for your analysis. This period should be long enough to accurately reflect the range of natural temperature variability. Commonly, a 30-year period is used in climate studies.

3. **Calculate Monthly Averages:** For each year in the baseline period, calculate the monthly average SST. This process involves averaging all daily temperature readings for each month.

4. **Identify the Maximum Monthly Average for Each Month:** For each month of the year, identify the highest monthly average SST recorded during the baseline period. This gives you 12 values, each representing the highest monthly average SST for a specific month across the baseline period.

5. **Determine the MMM:** The MMM is typically the highest of these 12 values. It represents the peak thermal threshold that the location experiences naturally.

6. **Consider Spatial Resolution:** Depending on your study's scale, you may need to consider the spatial resolution of your data. For large areas, you might need to calculate the MMM for different zones.

7. **Account for Climate Change and Variability:** Be aware that climate change can affect the MMM over time. Periodically updating the MMM with recent data can ensure that your analysis reflects current climate conditions.

8. **Use of Remote Sensing and Models:** In many cases, remote sensing data from satellites and oceanographic models are used to derive SST data, especially for remote or vast ocean areas.

9. **Verification and Validation:** Validate your calculated MMM with any available in-situ measurements or compare it with results from other studies or datasets to ensure accuracy.

By following these steps, you can determine the Historical Maximum Mean Temperature for a given location, which is crucial for assessing thermal stress on marine ecosystems, particularly coral reefs."

The [NOAA CoralReefWatch](https://coralreefwatch.noaa.gov/product/5km/tutorial/crw05a_sst_product.php) has tutorials on all their products.

I think the one I need to use is the [NCEI AVHRR Pathfinder SST](https://www.ncei.noaa.gov/products/avhrr-pathfinder-sst). The issue is, when I try to download the data files from 1985-2023, they are all .nc extensions. This is a NetCDF4 type file, which is not readable unless you use a code to convert it.

Some people in the field use MATLAB, but I'm going to try to use Python in the terminal.

ChatGPT: 

"Working with NetCDF4 data, which is a common format for storing large, multi-dimensional scientific data such as sea surface temperature (SST), can be done using various software tools. Here's a general guide on how to extract data from a NetCDF4 file into a more readable format:

### Tools and Libraries

1. **Python Libraries:** Python is a powerful tool for handling NetCDF4 data. Libraries like `netCDF4`, `xarray`, and `numpy` can be used to read, process, and visualize the data.

2. **GIS Software:** Geographic Information System (GIS) software like QGIS or ArcGIS can import NetCDF4 files directly and visualize the data on maps.

3. **MATLAB:** MATLAB is another popular tool for working with NetCDF4 data, especially in scientific and engineering fields.

### Steps for Python (Using `netCDF4` and `xarray`)

1. **Install Necessary Libraries:**
   ```python
   pip install netCDF4 xarray numpy matplotlib
   ```

2. **Import Libraries:**
   ```python
   import netCDF4 as nc
   import xarray as xr
   ```

3. **Open the NetCDF4 File:**
   ```python
   dataset = xr.open_dataset('path_to_your_file.nc')
   ```

4. **Explore the Data:**
   - List variables: `print(dataset.variables)`
   - Get info: `print(dataset.info())`

5. **Extract Data:**
   - Select specific data (e.g., a variable like SST): `sst_data = dataset['SST']` (replace 'SST' with the actual variable name)

6. **Convert to a Pandas DataFrame for Easier Reading (Optional):**
   ```python
   import pandas as pd
   df = sst_data.to_dataframe()
   ```

7. **Save to a CSV File (Optional):**
   ```python
   df.to_csv('output_filename.csv')
   ```

8. **Visualize the Data (Optional):**
   - Using `matplotlib` for plotting.

### Using GIS Software

1. **Open the GIS software and import the NetCDF4 file.**
2. **Navigate through the layers and select the appropriate SST data.**
3. **GIS tools allow you to visualize and even export the data into different formats.**

### MATLAB Approach

1. **Use MATLAB functions like `ncdisp`, `ncread`, or `netcdf.open` to explore and extract data from NetCDF4 files.**
2. **MATLAB also offers visualization tools to plot and analyze the data.**

### Additional Tips

- **Documentation:** Always refer to the documentation of the libraries or software you are using for specific functions and features.
- **Handling Large Data:** NetCDF4 files can be large; ensure your system has sufficient resources.
- **Learning Resources:** Online tutorials, forums, and community discussions can be helpful if you encounter specific issues or need advanced processing techniques.

By using these tools and steps, you can effectively convert the NetCDF4 data from NOAA's Coral Reef Watch into a more accessible and readable format for further analysis."

My code in the terminal:

```{bash}
python3 --version
pip3 install netCDF4 xarray numpy matplotlib #this line did not work because I need other packages (i.e. hdf5)

#install homebrew
 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

#install hdf5
brew install hdf5

#add this to ~/.bash_profile
export HDF5_DIR=$(brew --prefix hdf5)

export PATH=/Users/allysondemerlis/Library/Python/3.9/bin:$PATH #command line told me to do this

curl https://www.ncei.noaa.gov/data/oceans/pathfinder/Version5.3/

python3 #activate python

import netCDF4 as nc
import xarray as xr

dataset = xr.open_dataset('~/Downloads/19811231141339-NCEI-L3C_GHRSST-SSTskin-AVHRR_Pathfinder-PFV5.3_NOAA07_G_1981365_day-v02.0-fv01.0.nc')

sst_data = dataset['sea_surface_temperature'

sst_df = sst_data.to_dataframe()

print(sst_df.head())  # Displays the first 5 row

#longitude boundaries: -80.001, -80.17
#latitude boundaries: 25.60, 25.80

# Define your latitude and longitude range
lat_range = [26, 25]  # Replace with your desired range
lon_range = [-81, -80]  # Replace with your desired range

# Filter the data
filtered_data = sst_data.sel(lat=slice(*lat_range), lon=slice(*lon_range))

print(filtered_data)

import pandas as pd
df = filtered_data.to_dataframe()
df.to_csv('output_filename.csv')

```

**unresolved issue**: when I open this csv in R, it only contains data from 1981. WHY?
