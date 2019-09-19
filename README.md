# tpx3-to-root

Example software to demonstrate how to convert .tpx3 format to Root tree format.

## instructions

install and start root (https://root.cern.ch/)

.L tpx3_to_root.cpp+

tpx3_to_root("yourfile.tpx3")

this will generate yourfile.root

tpx3_clusters.cpp is used in the same way and will make a Root tree with cluster data.

The data folder contains two example .tpx3 files. One with data from a quad detector (512x512) and one from TPX3cam (256x256) data
with TDC events.

The variable kMaxPixel, found in tpx3_clusters.cpp, determines the maximum size of a stored cluster and can depend on the specific Timepix3 application.

## author
Erik Maddox (erik.maddoxATamscins.com), Amsterdam Scientific Instruments B.V.
