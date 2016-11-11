# Atmos-chem-models

Author: Jordan Capnerhurst 2016 B. Env. Science Honours student UOW

These scripts are my first attempt at coding so they're probably pretty sloppy and could be done a lot better but they may be of use to someone out there. Some are very straightforward and contain simple commands as I was learning python for the first time during the project.

The purpouse of them is to spatially and temporally plot NETCDF4 data outputted by 3 different chemical transport models.

The scripts included with this pacakage are for the assessment of the following chemical transport model outputs:

(KOEH)              CCAM-CTM 2011 
(Year OEH Scripts)  CCAM-CTM (V5.4) 
(KMEGAN)            CTM+MEGAN 2011 - Kathryn Emmerson CSIRO from paper "Current estimates of biogenic emissions from Eucalypts uncertain for Southeast Australia"
(JMEGAN/JMEGAN 3x3) MEGAN Offline 2011 - Jeremy Silver UOM (Both 1x1 and 3x3 data can be modelled from this)

These scripts assume that data has been concatanated into monthly files using CDO or other relevant software

Numerous scripts involve the datasets being plotted together eg. monthly timeseries. 
For these the data must be first imported using the individual model scripts and then combined by running the relevenat script in KMEGAN.

D1: llcrnrlon=147.804, llcrnrlat=-36.7246, urcrnrlon=153.114, urcrnrlat=-31.4146

D2: llcrnrlon=150.079, llcrnrlat=-34.836, urcrnrlon=151.849, urcrnrlat=-33.066

D3: llcrnrlon=150.1, llcrnrlat=-34.7177, urcrnrlon=151.651, urcrnrlat=-33.5651

(llcrn= lower left corner)
