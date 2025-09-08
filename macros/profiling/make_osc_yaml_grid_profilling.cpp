//  make_osc_yaml_grid.cpp  -----------------------------------------
//  Build :  g++ -std=c++17 -O2 make_osc_yaml_grid_profilling.cpp -lyaml-cpp -o make_osc_yaml_grid
//  Usage :  ./make_osc_yaml_grid         (writes 21 YAMLs for the parameter in `toScan`)
#include <yaml-cpp/yaml.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
namespace fs = std::filesystem;

// ------------------------------------------------------------------
const fs::path kTemplate = "../overrides/MaCh3_oscparameter_profiling.yaml";
const fs::path kOutDir   = "../outputs/profiling_yaml/MaCh3";
const std::string toScan = "PMNS_DELTA_MASS_SQUARED_32";

const std::unordered_map<std::string,double> kCentral = {
    {"PMNS_SIN_SQUARED_12",        0.310},
    {"PMNS_SIN_SQUARED_13",        0.0224},
    {"PMNS_SIN_SQUARED_23",        0.582},
    {"PMNS_DELTA_MASS_SQUARED_21", 7.39e-5},
    {"PMNS_DELTA_MASS_SQUARED_32", 0.002525},
    {"PMNS_DELTA_CP",              -2.498},
};

const std::vector<int> kSteps = {
    -10,-9,-8,-7,-6,-5,-4,-3,-2,-1,
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10};
// ------------------------------------------------------------------

// deep-clone helper
YAML::Node cloneDeep(const YAML::Node& n){
    YAML::Emitter e; e << n;
    return YAML::Load(e.c_str());
}

// return a copy of the parameter node
YAML::Node findParam(const YAML::Node& defs,const std::string& name){
    for(std::size_t i=0;i<defs.size();++i)
        if(defs[i]["parameterName"].as<std::string>()==name)
            return defs[i];
    throw std::runtime_error("Parameter "+name+" not found.");
}

int main(){
    if(!fs::exists(kTemplate)){
        std::cerr<<"Template "<<kTemplate<<" not found.\n";
        return 1;
    }

    YAML::Node tpl = YAML::LoadFile(kTemplate);
    fs::create_directories(kOutDir);

    // pull limits for the parameter we scan
    YAML::Node defsTpl =
        tpl["fitterEngineConfig"]["likelihoodInterfaceConfig"]
           ["propagatorConfig"]["parametersManagerConfig"]
           ["parameterSetList"][0]["parameterDefinitions"];
    YAML::Node nodeScan = findParam(defsTpl, toScan);

    double lo = nodeScan["parameterLimits"][0].as<double>();
    double hi = nodeScan["parameterLimits"][1].as<double>();
    double central = kCentral.at(toScan);
    double sigma   = (hi-lo)/10.1;

    std::cout<<"Profiling "<<toScan<<"  central="<<central
             <<"  sigma="<<sigma<<"\n";

    for(int n : kSteps){
        double newVal = central + n*sigma;
        if(newVal<lo || newVal>hi) continue;

        YAML::Node cfg  = cloneDeep(tpl);
        YAML::Node defs =
            cfg["fitterEngineConfig"]["likelihoodInterfaceConfig"]
               ["propagatorConfig"]["parametersManagerConfig"]
               ["parameterSetList"][0]["parameterDefinitions"];

        for(std::size_t i=0;i<defs.size();++i){
            std::string pname = defs[i]["parameterName"].as<std::string>();

            if(pname==toScan){
                defs[i]["isFixed"]          = true;
                defs[i]["priorValue"]       = newVal;
            }else{
                defs[i].remove("isFixed");
                defs[i]["priorValue"] = kCentral.at(pname);
            }
        }

        std::string tag   = (n>=0?"+":"") + std::to_string(n) + "sigma";
        fs::path   out    = kOutDir/("osc_"+toScan+"_"+tag+".yaml");

        std::ofstream fout(out);
        YAML::Emitter em; em.SetIndent(2); em << cfg;
        fout<<em.c_str();
        std::cout<<"→ "<<out<<"\n";
    }
    std::cout<<"Finished YAML generation in "
             <<fs::canonical(kOutDir)<<"\n";
    return 0;
}
