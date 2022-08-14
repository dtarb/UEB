# UEB

Utah Energy Balance Snowmelt Model

The Utah Energy Balance (UEB) snow model is an energy balance snowmelt model developed by David Tarboton's research group, first in 1994, and updated over the years. This repo is for the C++ version ported from the earlier Fortran version so as to be easier to use with NetCDF and parallel MPI libraries. The Fortran version is at [https://github.com/dtarb/UEBFortran](https://github.com/dtarb/UEBFortran)

UEB uses a lumped representation of the snowpack and keeps track of water and energy balance. The model is driven by inputs of air temperature, precipitation, wind speed, humidity and radiation at time steps sufficient to resolve the diurnal cycle (six hours or less). The model uses physically-based calculations of radiative, sensible, latent and advective heat exchanges. A force-restore approach is used to represent surface temperature, accounting for differences between snow surface temperature and average snowpack temperature without having to introduce additional state variables. Melt outflow is a function of the liquid fraction, using Darcy's law. This allows the model to account for continued outflow even when the energy balance is negative. Because of its parsimony (few state variables - but increasing with later versions) this model is suitable for application in a distributed fashion on a grid over a watershed. There are a number of versions available. 

### License

UEB is open source software: you can redistribute it and/or modify it under the terms of the MIT Open Source License as published by the Open Source Initiative https://opensource.org/licenses/MIT.

UEB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

If you wish to use or incorporate this program (or parts of it) into other software that does not meet the MIT License conditions contact the author to request permission.

David G. Tarboton  
Utah State University  
8200 Old Main Hill  
Logan, UT 84322-8200  
USA  
[http://hydrology.usu.edu/dtarb/](http://hydrology.usu.edu/dtarb/) 

email:  dtarb@usu.edu 

### Acknowledgements

I am grateful to the following funding agencies that have supported the develoopment of UEB

* NSF Grant EPS 1135482 CI-WATER for the development of the parallel C++ version
* NASA Grant NNX11AK03G has supported the development of UEBGrid
* USDA-CREES award 2008-34552-19042 Utah Drought Management Project supported the development of the vegetation components.
