# Example: Using NuOscillator to Tabulate Oscillations

This is a simple example of how to use a oscillation weight calculation with GUNDAM.  This example is based on NuOscillator https://github.com/dbarrow257/NuOscillator.git, but the basic principles should apply to any tabulation method.

Compiling this code on linux, you have to have a valid ROOT in your path which is typically done by sourcing the GUNDAM setup.sh, but can also be done using `thisroot.sh`. The `compile.sh` script will run those commands and is the best way to build:

```bash
./compile.sh
```

You can also compile by hand using:

```bash
mkdir build-$(uname -m)
cd build-$(uname -m)
cmake -DCMAKE_INSTALL_PREFIX=${PWD} ..
make install
```

That will grab the NuOscillator code, build it, and then build the
libTabulatedNuOscillator.so file that can be included in GUNDAM.  This
probably doesn't work on MacOS, but probably could be made to work.

> Note: Do to NuOscillator idiosyncracies, this needs to be installed into
> it's build directory.

After compiling, you need to source the NuOscillator setup script in `<build>/bin/setup.NuOscillator.sh`.

# Controlling TabulatedNuOscillator with the GUNDAM yaml config file

This library is used to implemented a `Tabulated` dial in the a gundam yaml
config file.  This is done under the each parameter set definition.  Here
is an example that setups the oscillation table for numu-to-numu,
numu-to-nue, and numubar-to-numubar.

```yaml
  fitterEngineConfig:

    propagatorConfig:

      parameterSetListConfig
      - name: "Oscillation Parameters"
        parameterDefinitions:
          - parameterName: "GUNDAM_SIN_SQUARED_12"
            isEnabled: true
            priorValue: 0.307
            priorType: Flat
            parameterStepSize: 1.0
            parameterLimits: [ 0.0, 1.0 ]
          - parameterName: "GUNDAM_SIN_SQUARED_13"
            isEnabled: true
            priorValue: 0.0219
            priorType: Flat
            parameterStepSize: 1.0
            parameterLimits: [ 0.0, 1.0 ]
          - parameterName: "GUNDAM_SIN_SQUARED_23"
            isEnabled: true
            priorValue: 0.558
            priorType: Flat
            parameterStepSize: 1.0
            parameterLimits: [ 0.0, 1.0 ]
          - parameterName: "GUNDAM_DELTA_MASS_SQUARED_21"
            isEnabled: true
            priorValue: 7.53E-5
            priorType: Flat
            parameterStepSize: 1.0
            parameterLimits: [ 0.0, 1.0 ]
          - parameterName: "GUNDAM_DELTA_MASS_SQUARED_32"
            isEnabled: true
            priorValue: 0.002455
            priorType: Flat
            parameterStepSize: 1.0
            parameterLimits: [ 0.0, 1.0 ]
          - parameterName: "GUNDAM_DELTA_CP"
            isEnabled: true
            priorValue: 3.74       # 1.19 pi
            priorType: Flat
            parameterStepSize: 1.0
            parameterLimits: [ 0.0, 2.0 ]

        dialSetDefinitions:

          - dialType: Tabulated
            applyCondition: "[GUNDAM_IS_CC]!=0&&[GUNDAM_FLUX]==14&&[GUNDAM_NU]==14"
            dialInputList:
              - name: "GUNDAM_SIN_SQUARED_12"
              - name: "GUNDAM_SIN_SQUARED_13"
              - name: "GUNDAM_SIN_SQUARED_23"
              - name: "GUNDAM_DELTA_MASS_SQUARED_12"
              - name: "GUNDAM_DELTA_MASS_SQUARED_32"
              - name: "GUNDAM_DELTA_CP_32"
            tableConfig:
              name: "Tabulated numu to numu NuFASTLinear"  # Must be unique.
              libraryPath: "${NUOSCILLATOR}/lib/libTabulatedNuOscillator.so"
              initFunction: "initializeTable"
              initArguments:
                - "CONFIG ${NUOSCILLATOR}/Configs/Unbinned_NuFASTLinear.yaml"
                - "FLUX_FLAVOR muon"          # muon neutrino flux
                - "INTERACTION_FLAVOR muon"   # muon neutrino
                - "BINNING_FILE exampleEnergyBinning.root"
                - "BINNING_HIST energyBinning"
                - "PATH 1300.0"
                - "PARAMETERS SS12,SS13,SS23,DM21,DM32,DCP"
              updateFunction: "updateTable"
              binningFunction: "bintable"
              binningVariables:
                - "GUNDAM_NU_ENERGY" # Nu energy in GeV
                - "GUNDAM_NU_ZENITH" # Nu zenith cosine (optional)

          - dialType: Tabulated
            applyCondition: "[GUNDAM_IS_CC]!=0&&[GUNDAM_FLUX]==-14&&[GUNDAM_NU]==-14"
            dialInputList:
              - name: "GUNDAM_SIN_SQUARED_12"
              - name: "GUNDAM_SIN_SQUARED_13"
              - name: "GUNDAM_SIN_SQUARED_23"
              - name: "GUNDAM_DELTA_MASS_SQUARED_12"
              - name: "GUNDAM_DELTA_MASS_SQUARED_32"
              - name: "GUNDAM_DELTA_CP_32"
            tableConfig:
              name: "Tabulated anti-numu to anti-numu NuFASTLinear"
              libraryPath: "${NUOSCILLATOR}/lib/libTabulatedNuOscillator.so"
              initFunction: "initializeTable"
              initArguments:
                - "CONFIG ${NUOSCILLATOR}/Configs/Unbinned_NuFASTLinear.yaml"
                - "FLUX_FLAVOR anti-muon"          # muon neutrino flux
                - "INTERACTION_FLAVOR anti-muon"   # muon neutrino
                - "BINNING_FILE exampleEnergyBinning.root"
                - "BINNING_HIST energyBinning"
                - "PATH 1300.0"
                - "PARAMETERS SS12,SS13,SS23,DM21,DM32,DCP"
              updateFunction: "updateTable"
              binningFunction: "bintable"
              binningVariables:
                - "GUNDAM_NU_ENERGY" # Nu energy in GeV
                - "GUNDAM_NU_ZENITH" # Nu zenith cosine (optional)

          - dialType: Tabulated
            applyCondition: "[GUNDAM_IS_CC]!=0&&[GUNDAM_FLUX]==14&&[GUNDAM_NU]==12"
            dialInputList:
              - name: "GUNDAM_SIN_SQUARED_12"
              - name: "GUNDAM_SIN_SQUARED_13"
              - name: "GUNDAM_SIN_SQUARED_23"
              - name: "GUNDAM_DELTA_MASS_SQUARED_12"
              - name: "GUNDAM_DELTA_MASS_SQUARED_32"
              - name: "GUNDAM_DELTA_CP_32"
            tableConfig:
              name: "Tabulated numu to nue NuFASTLinear"  # Must be unique.
              libraryPath: "${NUOSCILLATOR}/lib/libTabulatedNuOscillator.so"
              initFunction: "initializeTable"
              initArguments:
                - "CONFIG ${NUOSCILLATOR}/Configs/Unbinned_NuFASTLinear.yaml"
                - "FLUX_FLAVOR muon"              # muon neutrino
                - "INTERACTION_FLAVOR electron"   # electron neutrino
                - "BINNING_FILE exampleEnergyBinning.root"
                - "BINNING_HIST energyBinning"
                - "PATH 1300.0"
                - "PARAMETERS SS12,SS13,SS23,DM21,DM32,DCP"
              updateFunction: "updateTable"
              binningFunction: "bintable"
              binningVariables:
                - "GUNDAM_NU_ENERGY" # Nu energy in GeV
                - "GUNDAM_NU_ZENITH" # Nu zenith cosine (optional)
```

The initArguments values for everything except the `FLUX_FLAVOR`, and
`INTERACTION_FLAVOR` must be the same for all of the tables filled by
`libTabulatedNuOscillator.so`.  The neutrino flavor names are the strings
used by NuOscillator, with the prefix "anti-" specifying that the neutrino
type is an anti-neutrino.

| Neutrino Type         |        Flavor | PDG code |
|:----------------------|--------------:|---------:|
| electron neutrino     |      electron |       12 |
| electron antineutrino | anti-electron |      -12 |
| muon neutrino         |          muon |       14 |
| muon antineutrino     |     anti-muon |      -14 |
| tau neutrino          |           tau |       16 |
| tau antineutrino      |      anti-tau |      -16 |

The flux types are the same as the neutrino type.  When the flux type
matches the neutrino type, a survival probability is used.  When the flux
type does not match the neutrino type, an appearance probability is used.

The parameter definitions are contained in a comma separated list that
defines the order of the oscillation parameters being provided by GUNDAM.

## The parameter definitions

The parameter definitions are contained in a text file that specifies how
the GUNDAM config file has defined the parameters.  This is the order of
the parameters in the input array of doubles.  The oscillation parameters
are defined using the strings

* SS12    -- Sin-Squared theta 12.
* SS13    -- Sin-Squared theta 13.
* SS23    -- Sin-Squared theta 23.
* DM21    -- delta M^2 21.
* DM32    -- delta M^2 32.
* DCP     -- delta CP.

These values will be mapped into the enum values used by NuOscillator
(i.e. kTH12, kTH12, kTH23, kDM12, kDM23, and kDCP [Note: NuOscillator
reverses the indices on kDM23, so it is Delta M-squared 3 minus Delta
M-squared 2]).

## Oscillation weight grid

The oscillation weights will be calculated at values of neutrino energy and zenith angle (for atmospheric neutrino oscillations) that are provided using a histogram.  The histogram will be accessed as a TH1, TH2, or TH3 where the X-axis specifies the energy values to calculate, the Y-axis (if present) specifies the zenith cosine, and the Z-axis (if present) specifies the production height values [note: Most oscillators require that the production heights be uniformly spaced]. Example histograms can be constructed using the `exampleEnergyBinning.py`, `exampleZenithBinning.py`, and `exampleHeightBinning.py` scripts.  The grid is provided to the library using the `BINNING_FILE` and `BINNING_HIST` strings that are passed to the function specified by `initFunction`.

The grid values are accessed using the TAxis interface.  The values are either the low edge of the bin, the linear center of the bin `(H+L)/2`, the logarithmic center of the bin `exp((log(H)+log(L))/2)`, or the "inverse" center of the bin `(H-L)/(log(H)-log(L))` [this is the average `x` weighted by `1/x`, which is appropriate for "1/E" binning]. The NuOscillator binned oscillators use the linear center.  The way the grid is calculated from the histogram is set using `ENERGY_TYPE` (can be "edge", "average", "logarithmic", or "inverse"), and `ZENITH_TYPE` (can be "edge", or "average").

The contents of the TH1 and TH2 histograms used to specify the energy and zenith cosine distributions will be ignored.

When a TH3 is used to specify the production height distribution, the histogram contents will be the probability at each production height. In this case, the histogram will need to be normalized based on the needs of the oscillator (typically, that means the sum of each column must be normalized).  See `exampleHeightBinning.py` for a simple, non-physical, example that can be used with CUDAProb3.
