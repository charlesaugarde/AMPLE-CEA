# AMPLE.m
AMPLE - A Material Point Learning Environment

AMPLE.m is a quasi-static implicit implementation of the material point method in MATLAB.  
More information about AMPLE.m can be obtained from the project webapges:
https://wmcoombs.github.io/AMPLE

AMPLE is an elasto-plastic large deformation material point code with a regular quadrilateral background mesh 
(the main code is ample.m).   The continuum framework is based on an updated Lagrangian formation and two 
different constitutive models are included: linear elasticity and a linear elastic-perfectly plastic model 
with a von Mises yield surface.  
 
It is suggested that users start by looking at the “index.html” file within the “documentation” folder.  
The html files were complied with M2HTML (https://www.artefact.tk/software/matlab/m2html/) and they should 
provide  an overview of the code and how it all fits together.  The scripts should also be compatible with 
MATLAB’s “help” function. 

If you use AMPLE.m in your research, please cite the following paper in any associated publications:
Coombs, W., & Augarde, C. (2020). AMPLE: A Material Point Learning Environment. Advances in Engineering Software, 139, Article 102748. https://doi.org/10.1016/j.advengsoft.2019.102748

Version history:

AMPLE     - original AMPLE release

AMPLE_1.1 - version 1.1, released August 2020 - focused on runtime improvement

July 2026:
* added 3D elastic column collapse creating a new setup file called `setupGrid_column3D.m`
* for consistency renamed the 2D elastic column collapse to `setupGrid_column2D.m`
* created `formCoord3D.m` for the 3D example
* in testing the 3D problem came across convergence issues and generation of stresses in the two horizontal directions. Tracked this down to the generation of very small impaginary parts to `duvw`, caused by slight unsymmetry in `epsEN` and `BeT` not seen with the 2D version. Dealt with this by adding lines to ensure symmetry e.g. `epsEn = 0.5*(epsEn + epsEn');` in `detMPs.m`.
* to discuss with WMC - addition/amendment of commenting in the code, dates at the top etc.
  
August 2026
* added 1D column example `setupGrid_column1D.m`. Checked for agreement with 2D case for both standard MPM and GIMP
  
