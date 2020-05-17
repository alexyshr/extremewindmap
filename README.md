# extremewindmap
Spatio-temporal analysis of extreme wind velocities for infrastructure design. 
Dissertation submitted in partial fulfillment of the requirementsfor the Degree of 
Master of Science in Geospatial Technologies.

Alexys Herleym Rodríguez Avellaneda
alexyshr@gmail.com
https://github.com/alexyshr

Supervised by:
Prof. Dr. Edzer Pebesma
Institute for Geoinformatics
University of Münster - Germany

Co-supervised by:
Prof. Dr. Juan C. Reyes
Department of Civil and Environmental Engineering
Universidad de los Andes - Colombia

Co-supervised by:
Prof. Dr. Sara Ribero
Information Management School
Universidade Nova de Lisboa - Portugal

Abstract

This research aims to create non-hurricane non-tornadic maps of extreme wind speeds for
the mean recurrence intervals MRIs 700, 1700, and 3000 years, covering the Colombian
territory. For infrastructure design, these maps are combined with existing hurricane wind
speed studies, to be used as input loads due to wind.

For each station with non-thunderstorm wind speeds time histories in the input data, following
(Pintar, Simiu, Lombardo, & Levitan, 2015), extreme wind speeds corresponding to
each MRI are calculated using a Peaks Over Threshold Poisson Process POT-PP extreme
value model, then wind velocities with the same MRI are spatially interpolated to generate
continuous maps for the whole study area. The annual exceedance probability for all velocity
values in 700, 1700 and 3000 years MRIs output maps are respectively 1/700, 1/1700 and
1/3000.

Regarding input data, not only time series of field measurements from IDEAM methodological
stations are used, but also post-processed information coming from the Integrated Surface
Database ISD, and ERA5 forecast reanalysis data. This condition demanded a comparison
of the different data sources, in order to verify the feasibility in the use of ISD and ERA5,
this is downscaling support. The result of the comparison showed little similarity between
the different sources, but taking into account that complete and adequate measured data
from IDEAM was not available. Before to apply POT-PP, ISD and IDEAM data sources
were standardized to meet the requirement of three seconds (3-s) wind gust speed, ten (10)
meters anemometer height, and terrain open space condition.

Due to the limitation in the classification of thunderstorm and non-thunderstorm data, it was
not possible to take real advantage of POT-PP method, which was limited/restricted from
non-homogeneous to homogeneous and from non-stationary to stationary, being equivalent to
use the most common POT - generalized Pareto approach. Non-hurricane maps were created
for data sources ISD and ERA5, using Kriging as spatial interpolation method. After the
integration with previous hurricane studies, the results of ERA5 showed the most reliable
final maps, despite limitations in the input data to guarantee downscaling support. ISD final
map showed very high wind values, which are unlikely. These shortcomings may be corrected
when complete IDEAM data-source and storm data classification are available.

A complete R tool was implemented to solve the whole process, which is based in copyrighted
code for de-clustering and thresholding generously given by Dr Adam L. Pintar adam.pintar@nist.gov 
- National Institute of Standards and Technology NIST, U.S Department of Commerce.

Thesis document: https://github.com/alexyshr/extremewindmap/blob/master/_book/thesisoneside_final.pdf
