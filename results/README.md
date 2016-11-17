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