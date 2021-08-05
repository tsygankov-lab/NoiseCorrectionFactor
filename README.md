# NoiseCorrectionFactor


Dr. Denis Tsygankov (2021)
Integrative Systems Biology Lab, 
Wallace H. Coulter Department of Biomedical Engineering, Georgia Institute of Technology and Emory University SOM 
 
If you use any part of the scripts in this package, please cite:

DJ Marston et al. "Correcting artifacts in ratiometric biosensor imaging; an improved approach for dividing noisy signals"
Frontiers in Cell and Developmental Biology (2021), doi: 10.3389/fcell.2021.685825

The package includes two examples of data sets (from GEF Asef and CDC42 biosensors):

Numerator1.tif    (a multi frame signal for the numerator of the ratio)
Denominator1.tif (a multi frame signal for the numerator of the ratio)
Mask1.tif (a single frame image of the cell mask)

	and similarly 

Numerator2.tif    
Denominator2.tif 
Mask2.tif 

The scripts NCF_processing1.m and NCF_processing2.m illustrate the NCF method 
using the optimizations defined by Equations 5 and 6 in the paper, respectively.

If you use NCF_processing1.m for the first data set and NCF_processing2.m for the second data set,
you will see that results that correspond to Figures 8 and 9, respectively.

The package also includes profile.m, which was used to generate the profiles of the mean intensity along 
the cell perimeter as a function of the distance from the cell edge for Figures 1, 3, 4 and 6 in the paper.      
