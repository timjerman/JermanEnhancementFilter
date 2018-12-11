# Jerman Enhancement Filter
Jerman's 3D and 2D Hessian based tubular (vessel/vesselness) and spherical (blob/blobness) enhancement filters.

The MATLAB code is the implementation of the next two journal publications:

1. [T. Jerman, F. Pernus, B. Likar, Z. Spiclin, "*Enhancement of Vascular Structures in 3D and 2D Angiographic Images*", IEEE Transactions on Medical Imaging, 35(9), p. 2107-2118 (2016), doi={10.1109/TMI.2016.2550102}](https://doi.org/10.1109/TMI.2016.2550102)

2. [T. Jerman, F. Pernus, B. Likar, Z. Spiclin, "*Blob Enhancement and Visualization for Improved Intracranial Aneurysm Detection*", IEEE Transactions on Visualization and Computer Graphics, 22(6), p. 1705-1717 (2016), doi={10.1109/TVCG.2015.2446493}](https://doi.org/10.1109/TVCG.2015.2446493)

and a conference paper (please refer to journal publication [1] for more details):

3. [T. Jerman, F. Pernus, B. Likar, Z. Spiclin, "*Beyond Frangi: an improved multiscale vesselness filter*", Proc. SPIE 9413, Medical Imaging 2015: Image Processing, 94132A (2015), doi={10.1117/12.2081147}](https://doi.org/10.1117/12.2081147)

The code is based on Dirk-Jan Kroon's implementation of Frangi's vesselness filter. (https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter)

### Tips:

* Make sure that the objects of interest have the highest (if bright compared to the background) or lowest (if dark compared to background) intensities in the image/volume. Scale/normalize the images appropriately.

* The 3D method contains a c-code file that needs to be compiled with "mex eig3volume.c". (For more info visit: https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter)

* Threshold the filter response to remove any remaining enhanced noise

### Content:

2D enhancement of vessel/tube-like structures:

 * vesselness2D.m - main function
 * example_vesselness2D.m - filter applied on a 2D retinal vasculature
 * fundus2D.png - image for the example
 
3D enhancement of vessel/tube-like structures:

 * vesselness3D.m - main function
 * eig3volume.c - fast computation of eigenvalues
 * example_vesselness3D.m - filter applied on a 3D cerebral vasculature
 * volume.mat - volume for the example
 
3D enhancement of blob/sphere-like structures:

 * blobness3D.m - main function
 * eig3volume.c - (as above)
 * example_blobness3D.m - aneurysm enhancement in a 3D cerebral vasculature
 * volume.mat - (as above)
 
 2D enhancement of blob/disk-like structures (**experimental - not published**):

 * blobness2D.m - main function
 * example_blobness2D.m - enhancement of an aneurysm on a vessel (*synthetically generated image*)
 * blob2D.png - image for the example
  
