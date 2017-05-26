# Jerman Enhancement Filter
Jerman's 3D and 2D Hessian based tubular (vessel/vesselness) and spherical (blob) enhancement filters.

The MATLAB code is the implementation of the next two journal publications:

1. T. Jerman and F. Pernus and B. Likar and Z. Spiclin, "*Enhancement of Vascular Structures in 3D and 2D Angiographic Images*", IEEE Transactions on Medical Imaging, 35(9), p. 2107-2118 (2016), doi={10.1109/TMI.2016.2550102}

2. T. Jerman and F. Pernus and B. Likar and Z. Spiclin, "*Blob Enhancement and Visualization for Improved Intracranial Aneurysm Detection*", IEEE Transactions on Visualization and Computer Graphics, 22(6), p. 1705-1717 (2016), doi={10.1109/TVCG.2015.2446493}


The code is based on Dirk-Jan Kroon's implementation of Frangi's vesselness filter. (https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter)

The 3D method contains a c-code file that needs to be compiled with "mex eig3volume.c". (For more info visit: https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter)
