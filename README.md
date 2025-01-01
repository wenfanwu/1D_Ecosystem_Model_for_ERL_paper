# 1D_Ecosystem_Model_for_ERL_paper

This one-dimensional (1D) ecosystem model is developed primarily based on the [CoSiNE module](https://ccrm.vims.edu/schismweb/CoSiNE_manual_ZG_v5.pdf) of SCHISM modeling system. It is designed to study the response of ecosystem to different SPM conditions using a 1D ecosystem model. It is also part of Wu et al. (2023, Weak local upwelling may elevate the risks of harmful algal blooms and hypoxia in shallow waters during the warm season, *Environ. Res. Lett.* 18(11), 114031).

**CoSiNE_1D_main.m**: the core code of the ecosystem model.

**CoSiNE_1D_demo.m**: the demo used to demonstrate the usage of model.

**cosine.nml**: namelist file of the model configuration.

The model contains a total of 11 state variables, including two phytoplankton species (small phytoplankton and diatom), two zooplankton species (microzooplankton and mesozooplankton), four dissolved inorganic nutrients (nitrate, phosphate,ammonium, and silicate), two detritus organic matters (detritus nitrogen and detritus silicon), and dissolved oxygen. 

Author: Wenfan Wu, Postdoc, Virginia Institute of Marine Science, William & Mary

Email: [wwu@vims.edu](mailto:wwu@vims.edu)

