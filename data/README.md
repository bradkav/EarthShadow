### Data files

This folder contains a number of data files (cross sections, density profiles, etc.) which are loaded and used by the `EarthShadow` code. Further documentation will be added later on the precise contents of these files and how to make use of them outside of the code.


#### Differential Cross Sections

The folder `dcross_sections` contains *differential cross-sections* for dark matter-nucleus scattering tabulated in the files `dsigmadE_op_iso.dat`, where op= 1, 3, ..., 15 is one of the dark matter-nucleon effective operators we consider, and iso labels 9 isotopes in the Earth:

| iso        | Isotope Name      |
| ------------- |-------------  |
| 1 | Oxygen|
| 2 | Silicon|
|3 | Magnesium|
|4 | Iron |
|5 | Calcium|
|7| Sodium |
|8| Sulfur|
|9| Nickel|
|10| Aluminium|


Each file has 4 columns: 

|DM mass (GeV) | recoil energy (keV) | cs1 (cm^2 keV^-1) | cs2 (cm^2 keV^-1)|
| ------------- |-------------  |---|---|
| | | | |


The values cs1 and cs2 are such that for a given operator, isotope, dark matter mass, and recoil energy, the corresponding value of the differential dark matter-nucleus scattering cross-section is

dsigma/dEr = (1/v^2) (cs1+v_T^{\perp 2} cs2)

where the dark matter-nucleus relative velocity v is expressed in natural units, v_T^{\perp 2}=v^2-q^2/(4 mu^2) and mu is the dark matter-nucleus reduced mass.

In calculating the cross sections, the isovector coupling constants have been set to 0, and the isoscalar coupling constants to 2/mV^2, where mV=246.2 GeV is the electroweak scale. This corresponds to assuming the value 1/mV^2 for the coupling constants for protons and neutrons. Rescaling for other values of the coupling constants should then be straightforward.

#### Total Cross Sections

The folder `totcross_sections` contains *total cross-sections* for dark matter-nucleus scattering tabulated in the files `sigma_op_iso.dat`, where op labels the effective operator and iso labels the isotopes in the Earth (see above). 

Each file has 3 columns: 

|DM mass (GeV) | velocity (natural units)| cross-section (cm^2)|
|---|---|---|
| | | |

As in the case of the differential cross-sections, the isovector coupling constants have been set to 0, and the isoscalar coupling constants to 2/mV^2, where mV=246.2 GeV is the electroweak scale. This corresponds to assuming the value 1/mV^2 for the coupling constants for protons and neutrons. Rescaling for other values of the coupling constants should then be straightforward.
