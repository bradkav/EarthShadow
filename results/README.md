### Tables of results

#### Speed distributions

In the `HighMass` and `LowMass` folders, you will find a number of tabulated speed distributions, calculated taking into account the effects of Earth-scattering. Each file is named as `Speeddist_op=[operator]_mx=[A.BC]_gam=[X.YZ]_ps=0.10.dat`, where `[operator]` is either `1`, `8` or `12`. The results are for a scattering probability of 10%. The value `[X.YZ]` is the angle gamma in units of pi (so 1.00 is gamma = pi).

The results in the `HighMass` folder are for a 50 GeV DM particle (`[A.BC] = 50.00`) while the results in `LowMass` are for a 0.5 GeV DM particle (`[A.BC]=0.50`).

Each data file contains 4 columns:

| v [km/s]   | f(v) [s/km]    | f\_A(v) - f(v) [s/km]  | f\_D(v) [s/km] |
| --- | --- | --- | --- |

Here, f(v) is the free (unperturbed) speed distribution, f\_A(v) is the speed distribution after attenuation, and f\_D(v) is the enhancement due to deflection. The full, perturbed speed distribution is the sum of the last three columns.

#### Scattering probabilities

In the `pscat` folder, you will find tables giving the DM-nucleon couplings required for a given average scattering probability. The files are named as `couplings_op=[operator]_ps=[prob].dat`. As above, `[operator]` is either `1`, `8` or `12`. The scattering probability is then `[prob]` which can be `0.01`, `0.10` or `0.50`.

Each data file contains 2 columns:

| Log10[m\_x / GeV]   | c / GeV^(-2) |
| --- | --- |

Here, m\_x is the DM mass and c is the isoscalar DM-nucleon coupling for the operator in question.

#### Event Rates

In the `CRESST-II` folder are results for the ratio of the number of events in a CRESST-II-like detector with and without including the effects of Earth-Scattering. The files are named as `CRESSTRate_op=[operator]_mx=[A.BC]_ps=[prob].dat`. As above, `[operator]` is either `1`, `8` or `12`. So far, only the files with `[A.BC] = 0.50` and `[prob] = 0.10` are included. 

Each data file contains 4 columns:

| gamma/pi | N<sub>free</sub> - N<sub>A</sub> | N<sub>free</sub> + N<sub>D</sub> | N<sub>free</sub> - N<sub>A</sub> + N<sub>D</sub> |
| --- | --- | --- | --- |

Here, gamma is the angle of the detector with respect to the DM wind, N<sub>free</sub> is the number of events without Earth-Scattering, N<sub>A</sub> is the contribution from attenuation and N<sub>D</sub> the contribution from deflection. The results are normalised such that N<sub>free</sub> = 1.