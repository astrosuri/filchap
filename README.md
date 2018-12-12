# Filament Characterization Package
FILament CHAracterization Package (FilChaP) is a python-based algorithm that can be employed to derive filament properties.
FilChaP is publicly available through this GitHub page. 
If you use FilChaP in a publication, please cite Suri et al. 2018.

https://doi.org/10.5281/zenodo.2222325


## Dependencies 

FilChaP is compatible with python 2.7 and requires the following libraries:

* Numpy,
* Astropy,
* Matplotlib,
* Pandas,
* Scipy.

We are currently working on making the code compatible with Python 3, please look out for new releases!

## Input & Output files

FilChaP takes in 2D (position-position) or 3D datasets (position-position-velocity) in FITS format,
and the coordinates of the filaments identified on the data as a text file. The organization of the text file can be found
in the example folder. The example text file is an output of DisPerSE (Sousbie 2011). In the example folder, I also a provide
a filament extraction algorithm. In case the filaments are identified with DisPerSE and extracted as a FITS file, one can use 
the provided algorithm here to convert the information in this FITS file to individual coordinate files for each filament. 
The individual coordinate files are then named 'filamentX.txt', where X is the filament ID. Th filament coordinate files include
x, y, z (if in PPV) coordinates of each point along the filament. FilChaP needs these coordinates to calculate the filament 
properties.

FilChaP results are large numpy arrays. For example, the **calculateWidth** function (see below) returns 32 parameters for each
slice along the filament. If there are 100 slices along a filament, the result array in this case has a size of 32x100.
The results arrays can be written in text files, but in order to deal with such large arrays, we recommend turning on the
pandas option, so that the results are saved in a HDF file in binary format. HDF files are also easier to query or plot.

## User guide

The following are the following functions that can be called within FilChaP:

* **calculateWidth** - A function to calculate filament widths. 
* **calculateLength** - A function to calculate the length of a filament.
* **calculateCurvature** - A function to calculate filament curvature.

### User defined parameters

FilChaP has two different set of parameters: (1) those that can be tuned by exploring the parameter space, (2) those that are the
properties of the data set and should be absolute values.

#### First category

These are the parameters that can be explored by the user:

* **npix** - the number of pixels that will be used for perpendicular cuts. The resulting
length of the perpendicular cut is equal to twice the npix value. If nothing is
known about the width of the filaments prior to running FilChaP , one can give a
rough estimation of about a few times 0.1 pc.
* **avg_len** - the length (in pc) along the filament over which the intensity profiles should be averaged. 
A meaningful number is 3 times the beamsize of the observations.
* **niter** - number of iterations for the baseline subtraction process of the intensity.
profiles
* **lam** - weight for baseline subtraction.
* **smooth** - number of pixels to smooth the intensity profile with. This is needed to
calculate the positions of the intensity maxima and minima which are needed to
derive the positions of fitting boundaries.

#### Second category

These are the parameters that are usually specific to the dataset.

* **fits_file** - path to the fits file that will be used to derive intensities or column
densities.
* **distance** - distance to the source in units of pc. Needed to calculate lengths within
the cloud.
* **pixel_size:** - size of 1 pixel in units of arcsec. Also needed to calculate lengths
within the cloud.
* **dv** - velocity resolution of the data in km s−1 , if it is in PPV. This is needed when
averaging the intensity profiles. 
* **noise_level** - the noise level of the average intensity profile to calculate reduced
chi-sq value. I note that, although the noise level is specific to an individual data
set, it varies with the changing number of averaged profiles. Therefore, the user
has to adjust this parameter whenever the avg_len is changed.


## Contact

Sümeyye Suri <br /> 
Max Planck Institute for Astronomy, <br /> 
Königstuhl 17, <br /> 
D-69117 Heidelberg, <br /> 
Germany

suri at mpia.de
