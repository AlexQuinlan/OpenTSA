CFRP CYCLIC DATA

Matlab files (.mat) - structures
Use Matlab function 'load'.

LF=[30.1 20.1 10.1 5.1 3.1] - loading frequency
Specimens: 
- CP0 - [0,90]_{3s}
- CP90 - [90,0_{3s}
- S45 - [+/-45]_{3s}
TSA data is divided in two datasets:
1) for mean values across the specimen surface (CFRP_TSA_mean.mat) - vector:
	For both cross-ply samples: CFRP.(sample).FR383.(parameter) where sample is CP0 or CP90; and parameter is Tmean, dTmean, pmean, dTstd, Tstd, pstd for mean temperature, mean peak-to-peak temperature, mean phase, standard deviation of peak-to-peak temperature, stardard deviation of mean temperature, standard deviation of phase. - mean values across sample surface.
	For [+/-45]_{3s} sample: CFRP.S45.FR383.amp900.(parameter) where parameter is the same as above.

2) full-field data of the specimen surface (CFRP_TSA_fields.mat) - array:
	CFRP.(sample).FR383.(LF).(parameter) where sample is CP0, CP90 or S45; LF is LF301, LF201, LF101, LF51 or LF31 (corresponding to 30.1, 20.1, 10.1, 5.1, and 3.1 Hz respectively); and parameter is T0 (mean temperature), dT (peak-to-peak temperature) or p (phase)


The LSDIC dataset includes both mean and full-field data (after the DIC processing with DaVis software and least-squares to extract the peak-to-peak strains)
LSDIC Exx correspond to the strain transverse to the load and Eyy to the strain longitudinal to the load (interchanged with respect to the nomenclature on the paper! Exx is longitudinal to the load and Eyy is transverse to the load)
(LSDIC_CFRP_mean.mat)
1) For full-field data  - array:
	LSDIC.(sample).(LF).(strain).(parameter) where sample is CP0, CP90 or S45; strain is Exx, Eyy or Exy in strain; LF is LF301, LF201, LF101, LF51 or LF31 (corresponding to 30.1, 20.1, 10.1, 5.1, and 3.1 Hz respectively); parameter is FieldAmp for peak-to-peak strain, FieldPhase for phase value of the strain with respect the first image or FieldOffset for mean strain (the one used in the paper is FieldAmp).
2) For mean value across the sample surface - vector:
	LSDIC.(sample).(strain).vector.(parameter) where parameter is meanval (mean value across the sample surface of the strain selected) or stdval (standard deviation across the sample surface of the strain selected) sorted for LF from 30.1 to 3.1 Hz

Strain gauges dataset includes the peak-to-peak strain extracted with the least-squares algorithm - vector:
(StrainGauges_CFRP.mat)
	SG.(sample).(parameter) where parameter is Exx or Eyy for peak-to-peak values of Exx and Eyy (being Exx longitudinal and Eyy transverse to the load - as in paper) in strains; or Exxstd and Eyystd for standard deviation 
	The strains are sorted for loading frequencies from 30.1 Hz to 3.1 Hz.

RAW DATA folder includes test machine and strain gauges data in .xlsx files - there is another read me file inside that folder