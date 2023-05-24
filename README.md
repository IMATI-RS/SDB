# SDB 

## A Matlab script to compute Satellite Derived Bathymetry (SDB) from Copernicus Satellite S-2

The Script_SDB.m is a MATLAB code to calculate Satellite Derived Bathymetry accordin to the band ratio technique improved by Stumpf et al. (Stumpf, R.P.; Holderied, K.; Sinclair, M. Determination of water depth with high-resolution satellite imagery over variable bottom types. Limnology and Oceanography 2003, 48, 547–556) 

The script requires as input the preprocessed satellite image product and the set of points for calibration and validation. 

In the following, a usage guide for the script is provided, divided into A) pre-processing of the satellite image (Copernicus Sentinel-2a) using the open SNAP tool, which involves the use of a batch graph for the relevant steps, which generates a pre-processed image, which is input to the B) processing that produces the derived bathymetry as .tif and .png files, along with some graphs and .csv files containing statistics. 

The guide takes as an example a use case related to Venice Lagoon, whose reference data are contained in the Data folder, while the Sentinel-2 product of the area is downloadable from ESA's Open Access Hub at the link (https://scihub.copernicus.eu/dhus/odata/v1/Products('980473c2-f2a3-4926-8827-9f57e8f919af')/$value).

As a final note, we note that the pre-processing step can be bypassed by providing as input to the script a georeferenced image with UTM/WGS84 projection, with the 4 bands Blue, Green, Red, and NIR.

## Requirements
- SNAP
- MatLab 


## A) Pre-processing (SNAP)
The early steps of the optical satellite image pre-processing procedure are contained in an .xml file that should be uploaded and executed on the SNAP desktop app.

NOTE: Before opening the products user should account to disable the read of aux files since some problems are related in Deglint phase (4) to these data: Deselect Sentinel-2 auxiliary data read, Apply, OK.

![Image](Figures/AUX_DATA_READ.png?raw=true)

1.	Drag and drop the optical satellite product zip file (S2A_MSIL2A_20200319T101021_N0214_R022_T32TQR_20200319T130518) in the SNAP product window 
2.	Open RGB image with the right click on the product for a visual check [B4, Red, B3, Green, B2, Blue for Sentinel-2 products]
3.	Open the graph builder tool from the SNAP toolbar and load the .xlm Graph file (Preprocessing_Graph), the following setting should be perfomed:

![Image](Figures/GRAPH.png?raw=true)

  - Select the desired input product from the available list ([1] S2A_MSIL2A_20200319T101021_N0214_R022_T32TQR_20200319T130518)
  - Resample the product by “reference band” and choose “B2” (or other 10 m resolution bands)
  - Subset according to:
    - The source bands to be used in the following steps (B2, B3, B4, B8 to select);
    - The desired AOI polygon, leveraging on Wicket -  https://arthur-e.github.io/Wicket/sandbox-gmaps3.html for the generation of the geographic coordinates in a proper format to be pasted in the proper space

The following steps are meant to:
  - Calculate the Land mask 
  - Merge the produced bands (User should select the five bands available in the window)
  - Apply the Land mask to the four source bands
  - Merge the produced bands (User should select the four bands available in the window)

![Image](Figures/DEGLINT.png?raw=true)

  - Save the pre-processed product (S2A_MSIL2A_20200319T101021_N0214_R022_T32TQR_20200319T130518_resampled_Subset_BandMaths) in BEAM-DIMAP format in the Data folder

4.	Deglint 
  - Create glint polygons or load polygons (Glint_Polygon) as ESRI shapefile: Vectos > Import > ESRI Shapefile; import the polygons as unique mask.
  - Open the Deglin processor: Optical > Thematic water processing > Sen2Coral > Deglint Processor   
  - Read the pre-processed product in the “Input/Output” parameters window ([2] S2A_MSIL2A_20200319T101021_N0214_R022_T32TQR_20200319T130518_resampled_Subset_BandMaths)
  - Save as GeoTIFF-BigTIFF
  - In the “processing parameters” write the polygons shapefile name for the “Sun Glint Areas”, select Blue, Green and Red land-masked bands as source, NIR as reference
  - Check the boxes
  - Run the processor
  - Read the product
  - Read vector as single layer
  - Write the polygons name
  - Select the bands to be corrected (B2, B3, B4) and the reference bands (B8, B9)
  - Check box for the
  - Write the product in the Preprocessed folder selecting GeoTiff/BigTiff format;

The batch processor can be used with the same procedure for a multiple images processing or in case of random errors with the graph builder
5.	End
It is also possible to perform the individual steps described in the same order by calling up each action from the SNAP toolbar.


 
## B) BATHRYMETRY DERIVATION (MATLAB)
Instructions for starting the code are included in the script. The work folder should contain the code file together with the shapefile of the ground truth points and the pre-processed satellite image.
 
## APPENDIX - ACCESS TO SENTINEL-2 SATELLITE DATA
Free registration on the Copernicus Data Space Ecosystem to access and download Sentinels’ collections - https://identity.dataspace.copernicus.eu/auth/realms/CDSE/login-actions/registration?client_id=cdse-public&tab_id=zXE8Re-jUkQ

Instructions and general information- https://www.youtube.com/watch?v=Am93Xi0PZ5o&ab_channel=CopernicusEU

Copernicus Ecosystem Documentation - https://documentation.dataspace.copernicus.eu/

Once registered and then in possession of credentials, the user accesses the SEARCH menu of the Data Browser to select:

1. Data source (e.g. Sentinel-2 - L2A)
2. Time Range
3. Cloud cover accepted percentage
4. The area of interest (AOI) by drawing a polygon on the online viewer with the appropriate tool.	

User can then download the desired product from the list.

 
 

