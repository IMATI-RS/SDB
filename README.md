# SDB 
##A Matlab script to compute Satellite Derived Bathymetry (SDB) from Copernicus Satellite S-2

SATELLITE DERIVED BATHYMETRY ALGORITHM 
IMAGE PREPROCESSING – SNAP 
The early steps of the optical satellite image pre-processing procedure are contained in an .xml file that should be uploaded and executed on the SNAP desktop app.
NOTE: Before opening the products user should account to disable the read of aux files since some problems are related in Deglint phase (4) to these data: Deselect Sentinel-2 auxiliary data read, Apply, OK.

![Image](Figures/AUX_DATA_READ.png?raw=true)

1.	Drag and drop the optical satellite product zip file (S2A_MSIL2A_20200319T101021_N0214_R022_T32TQR_20200319T130518) in the SNAP product window 
2.	Open RGB image with the right click on the product for a visual check [B4, Red, B3, Green, B2, Blue for Sentinel-2 products]
3.	Open the graph builder tool from the SNAP toolbar and load the .xlm Graph file (Preprocessing_Graph), the following setting should be perfomed:

![Image](Figures/GRAPH.png?raw=true)

a.	Select the desired input product from the available list ([1] S2A_MSIL2A_20200319T101021_N0214_R022_T32TQR_20200319T130518)
b.	Resample the product by “reference band” and choose “B2” (or other 10 m resolution bands)
c.	Subset according to:
i.	The source bands to be used in the following steps (B2, B3, B4, B8 to select);
ii.	The desired AOI polygon, leveraging on Wicket -  https://arthur-e.github.io/Wicket/sandbox-gmaps3.html for the generation of the geographic coordinates in a proper format to be pasted in the proper space
The following steps are meant to:
d.	Calculate the Land mask 
e.	Merge the produced bands (User should select the five bands available in the window)
f.	Apply the Land mask to the four source bands
g.	Merge the produced bands (User should select the four bands available in the window)

![Image](Figures/DEGLINT.png?raw=true)


h.	Save the pre-processed product (S2A_MSIL2A_20200319T101021_N0214_R022_T32TQR_20200319T130518_resampled_Subset_BandMaths) in BEAM-DIMAP format in the Data folder

4.	Deglint 
a.	Create glint polygons or load polygons (Glint_Polygon) as ESRI shapefile: Vectos > Import > ESRI Shapefile; import the polygons as unique mask.
b.	Open the Deglin processor: Optical > Thematic water processing > Sen2Coral > Deglint Processor   
c.	Read the pre-processed product in the “Input/Output” parameters window ([2] S2A_MSIL2A_20200319T101021_N0214_R022_T32TQR_20200319T130518_resampled_Subset_BandMaths)
d.	Save as GeoTIFF-BigTIFF
e.	In the “processing parameters” write the polygons shapefile name for the “Sun Glint Areas”, select Blue, Green and Red land-masked bands as source, NIR as reference
f.	Check the boxes
g.	Run the processor

h.	Read the product
i.	Read vector as single layer
j.	Write the polygons name
k.	Select the bands to be corrected (B2, B3, B4) and the reference bands (B8, B9)
l.	Check box for the
m.	Write the product in the Preprocessed folder selecting GeoTiff/BigTiff format;

The batch processor can be used with the same procedure for a multiple images processing or in case of random errors with the graph builder
5.	End
It is also possible to perform the individual steps described in the same order by calling up each action from the SNAP toolbar.


 
BATRHYMETRY DERIVATION - MATLAB
Instructions for starting the code are included in the script. The work folder should contain the code file together with the shapefile of the ground truth points and the pre-processed satellite image.
 
APPENDIX - ACCESS TO SENTINEL-2 SATELLITE DATA
Free registration on the Copernicus Data Space Ecosystem to access and download Sentinels’ collections - https://identity.dataspace.copernicus.eu/auth/realms/CDSE/login-actions/registration?client_id=cdse-public&tab_id=zXE8Re-jUkQ
Instructions and general information- https://www.youtube.com/watch?v=Am93Xi0PZ5o&ab_channel=CopernicusEU
Copernicus Ecosystem Documentation - https://documentation.dataspace.copernicus.eu/
Once registered and then in possession of credentials, the user accesses the SEARCH menu of the Data Browser to select:
1. Data source (e.g. Sentinel-2 - L2A)
2. Time Range
3. Cloud cover accepted percentage
4. The area of interest (AOI) by drawing a polygon on the online viewer with the appropriate tool.	
User can then download the desired product from the list.

 
 

