#include <iostream>
#include <stdexcept>
#include <string>

#include <yaml-cpp/yaml.h>

#include "TFile.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TObjArray.h"
#include "TObjString.h"

int main(int argc, char** argv) 
{
    const std::string yamlFile = argv[1];

    const int N = 208;

    YAML::Node root = YAML::LoadFile(yamlFile);
    YAML::Node systematicsHdr = root["Systematics"];

    std::string parnameStr = "xsec_param_names";
    TObjArray* parnameArray = new TObjArray();
    parnameArray->SetOwner(kTRUE);
    parnameArray->SetName(parnameStr.c_str());

    std::string nupdgStr = "xsec_norm_nupdg";
    TObjArray* nupdgArray = new TObjArray();
    nupdgArray->SetOwner(kTRUE);
    nupdgArray->SetName(nupdgStr.c_str());

    TVectorT<Double_t> priorvalVec(N);
    TVectorT<Double_t> nomvalVec(N);
    TVectorT<Double_t> paramvalLbVec(N);
    TVectorT<Double_t> paramvalUbVec(N);
    TVectorT<Double_t> kinematicLbVec(N);
    TVectorT<Double_t> kinematicUbVec(N);
    TVectorT<Double_t> stepscaleVec(N);

    TVectorD* nupdgVec = new TVectorD(N);

    TMatrixT<Double_t> covmat(N, N);
    std::vector<double> errorVec;

    for (int i = 0; i < N; ++i) 
    {
        double errval = systematicsHdr[i]["Systematic"]["Error"].as<double>();
        errorVec.push_back(errval);
    }   

    for (int i = 0; i < N; ++i)
    {
        YAML::Node systematic = systematicsHdr[i]["Systematic"];

        std::string paramName = systematic["Names"]["ParameterName"].as<std::string>();
        TObjString* s = new TObjString();
        s->SetString(paramName.c_str());
        parnameArray->Add(s);

        (*nupdgVec)[i] = systematic["NeutrinoFlavourUnosc"][0].as<double>();

        priorvalVec[i] = systematic["ParameterValues"]["PreFitValue"].as<double>();
        nomvalVec[i] = systematic["ParameterValues"]["Generated"].as<double>();
        paramvalLbVec[i] = systematic["ParameterBounds"][0].as<double>();
        paramvalUbVec[i] = systematic["ParameterBounds"][1].as<double>();
        kinematicLbVec[i] = systematic["KinematicCuts"][0]["TrueNeutrinoEnergy"][0].as<double>();
        kinematicUbVec[i] = systematic["KinematicCuts"][0]["TrueNeutrinoEnergy"][1].as<double>();
        stepscaleVec[i] = systematic["StepScale"]["MCMC"].as<double>();

        YAML::Node corrList = systematic["Correlations"];

        for (int j = 0; j < N; ++j) 
        {
            YAML::Node entry = corrList[j];
            double corrval = 0.0;
            for (auto it = entry.begin(); it != entry.end(); ++it) 
            {
                corrval = it->second.as<double>();
                break;
            }

            covmat(i, j) = corrval*errorVec[i]*errorVec[j];
        }
    }

    nupdgArray->Add(nupdgVec);

    // --------------------------------- Storage in ROOT file ---------------------------------//

    TFile outFile("flux_covariance_DUNE_ND+FD.root", "RECREATE");

    outFile.WriteObject(parnameArray, parnameStr.c_str());
    outFile.WriteObject(nupdgArray, nupdgStr.c_str());

    priorvalVec.Write("xsec_param_prior");
    nomvalVec.Write("xsec_param_nom");
    paramvalLbVec.Write("xsec_param_lb");
    paramvalUbVec.Write("xsec_param_ub");
    kinematicLbVec.Write("xsec_norm_kinematic_lb");
    kinematicUbVec.Write("xsec_norm_kinematic_ub");
    stepscaleVec.Write("xsec_stepscale");  

    covmat.Write("xsec_cov");

    outFile.Close();

    return 0;
}

