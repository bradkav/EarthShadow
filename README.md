# EarthShadow

### Tables of results

In the 'results' folder, you will find a number of tabulated results. Each file is named as `Speeddist_op=[operator]_mx=0.50_gam=[X.YZ]_ps=0.10.dat`, where `[operator]` is either `1`, `8` or `12`. The results are for a DM mass of 0.5 GeV and a scattering probability of 10%. The value `[X.YZ]` is the angle gamma in units of pi (so 1.00 is gamma = pi).

Each data file contains 4 columns:

| v [km/s]   | f(v) [s/km]    | f\_A(v) - f(v) [s/km]  | f\_D(v) [s/km] |
| --- | --- | --- | --- |

Here, f(v) is the free (unperturbed) speed distribution, f\_A(v) is the speed distribution after attenuation, and f\_D(v) is the enhancement due to deflection. The full, perturbed speed distribution is the sum of the sum of the last three columns.
