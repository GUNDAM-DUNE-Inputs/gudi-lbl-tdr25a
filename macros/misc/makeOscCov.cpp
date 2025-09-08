//  makeOscCov.cpp  -----------------------------------------------------------
//  Run, e.g.:
//      root -l -b -q 'makeOscCov.cpp("./inputs/parameters/osc_cov.root")'

#include <TFile.h>
#include <TMath.h> 
#include <TObjArray.h>
#include <TObjString.h>
#include <TMatrixT.h>
#include <iostream>

void makeOscCov(const char* outfile = "./inputs/parameters/osc_cov.root")
{
    // ------------------------------------------------------------------ //
    // 1. Parameter names
    // ------------------------------------------------------------------ //
    const char* names[] = {
        "PMNS_SIN_SQUARED_12",
        "PMNS_SIN_SQUARED_13",
        "PMNS_SIN_SQUARED_23",
        "PMNS_DELTA_MASS_SQUARED_21",
        "PMNS_DELTA_MASS_SQUARED_32",
        "PMNS_DELTA_CP"
    };
    constexpr int N = sizeof(names) / sizeof(names[0]);

    auto* nameArray = new TObjArray(N);
    nameArray->SetName("osc_param_names");
    for (int i = 0; i < N; ++i)
        nameArray->Add(new TObjString(names[i]));

    // ------------------------------------------------------------------ //
    // 2. 6×6 diagonal covariance matrix
    // ------------------------------------------------------------------ //
    // const double diagVals[N] = {0.013, 0.0007, 0.021, 1.8e-6, 3.4e-5, 1.17}; // MaCh3
    const double diagVals[N] = {0.013, 0.0007, 0.021, 1.8e-6, 2.8e-5, 0.691}; // PDG: https://pdg.lbl.gov/2024/listings/rpp2024-list-neutrino-mixing.pdf :Three-neutrino mixing parameters
    

    auto* cov = new TMatrixT<double>(N, N);
    cov->Zero();
    for (int i = 0; i < N; ++i){
        // (*cov)(i, i) = TMath::Sqrt(diagVals[i]);
        (*cov)(i, i) = diagVals[i]*diagVals[i];
    }
        
    // ------------------------------------------------------------------ //
    // 3. Write to file
    // ------------------------------------------------------------------ //
    TFile f(outfile, "RECREATE");

    //  Write the array as ONE key only:
    //     kSingleKey = do NOT write the six sub-objects as individual keys
    nameArray->Write("osc_param_names", TObject::kSingleKey);

    // The matrix still goes under its own key
    f.WriteObject(cov, "osc_cov");

    f.Close();

    std::cout << "Wrote " << outfile << " with:\n"
              << "  • osc_param_names (TObjArray, single key)\n"
              << "  • osc_cov         (TMatrixT<double>, 6×6)\n";
}
