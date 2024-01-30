# NEMO4.2_P6Z

PISCES QUOTA P6Z includes explicit diazotrophy for the standard P5Z PISCES QUOTA framework

Able to switch between two prevalent marine diazotrophs, Trichodesmium and Crocosphaera

Able to account for the temperature dependence of the nitrogen fixation elemental use efficiencies of both Iron (Fe) and Phosphorus (P) for both diazotrophs

Nitrogen fixation is facultative within the model

Diazotroph thermal performance curves for growth and elemental use efficiencies based on Jiang et al. (2018) and Yang et al. (2021) for Trichodesmium and Crocosphaera, respectively

Jiang, H.-B., Fu, F.-X., Rivero-Calle, S., Levine, N. M., Sañudo-Wilhelmy, S. A., Qu, P.-P., Wang, X.-W., Pinedo-Gonzalez, P., Zhu, Z., & Hutchins, D. A. (2018). Ocean warming alleviates iron limitation of marine nitrogen fixation. Nature climate change, 8(8), 709-712. doi:10.1038/s41558-018-0216-8

Yang, N., Merkel, C. A., Lin, Y.-A., Levine, N. M., Hawco, N. J., Jiang, H.-B., Qu, P.-P., DeMers, M. A., Webb, E. A., Fu, F.-X., & Hutchins, D. A. (2021). Warming iron-limited oceans enhance nitrogen fixation and drive biogeographic specialization of the globally important cyanobacterium Crocosphaera. Frontiers in Marine Science, 8(118). doi:10.3389/fmars.2021.628363

MY_SRC (Model subroutines)

EXPREF directory (namelists, file_def, field_def and other xml files)

_In namelist_pisces_ref_p6z:_

Select flag ln_p6z to select explicit diazotrophy

Flag ln_tricho chooses diazotroph setup: _Trichodesmium_ (True) and _Crocosphaera_ (False)

Flag ln_facul activates facultative diazotrophy 

Flag ln_tiue controls the temperature dependence of the nitrogen fixation Fe use efficiency 

Flag ln_tiue controls the temperature dependence of the nitrogen fixation P use efficiency 

REFERENCE CONFIGURATION:

ln_tricho = true
ln_facul = true
ln_tiue = false
ln_piue = false
