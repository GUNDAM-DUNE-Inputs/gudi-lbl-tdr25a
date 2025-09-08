// profiling_scan_plot.C  --------------------------------------------------
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <filesystem>
#include <regex>
#include <map>
#include <vector>
#include <algorithm>
namespace fs = std::filesystem;

const char*  TREEPATH  = "FitterEngine/postFit/bestFitStats";
const char*  LEAF_LLH  = "totalLikelihoodAtBestFit";
const fs::path PLOT_DIR = "outputs/profiling_plots/MaCh3";

/* filename pattern: osc_<PARAM>_+3sigma  or  osc_<PARAM>__3sigma */
const std::regex FILE_RE("osc_([^_]+(?:_[^_]+)*)_([+]?|__?)(\\d)sigma");

/* parameter → leaf index and pretty axis label                   */
struct LeafInfo { int idx; const char* label; };
const std::map<std::string,LeafInfo> leafInfo = {
  // {"PMNS_SIN_SQUARED_12",        {0, "sin^{2}#theta_{12}"}},
  // {"PMNS_SIN_SQUARED_13",        {1, "sin^{2}#theta_{13}"}},
  // {"PMNS_SIN_SQUARED_23",        {2, "sin^{2}#theta_{23}"}},
  // {"PMNS_DELTA_MASS_SQUARED_21", {3, "#Deltam^{2}_{21} (eV^{2})"}},
  {"PMNS_DELTA_MASS_SQUARED_32", {4, "#Deltam^{2}_{32} (eV^{2})"}},
  // {"PMNS_DELTA_CP",              {5, "#delta_{CP} (rad)"}}
};

void profiling_scan_plot(const char* dir = "outputs/profiling_root/MaCh3")
{
  if (!fs::exists(dir)) { std::cerr<<"Dir "<<dir<<" not found\n"; return; }
  fs::create_directories(PLOT_DIR);

  std::map<std::string,std::vector<std::pair<double,double>>> pts;  // par→(x,y)

  for (const auto& ent : fs::directory_iterator(dir))
  {
    if (ent.path().extension() != ".root") continue;

    std::smatch m;
    std::string stem = ent.path().stem().string();
    if (!std::regex_search(stem, m, FILE_RE)) continue;

    std::string par = m[1];
    auto lf = leafInfo.find(par);
    if (lf == leafInfo.end()) continue;                   // unknown param

    TFile f(ent.path().c_str(),"READ");
    auto* t = static_cast<TTree*>(f.Get(TREEPATH));
    if (!t) continue;

    // ---- y value (likelihood)
    double llh = 0.0;

    t->GetEntry(0);
    llh = t->GetLeaf(LEAF_LLH)->GetValue();

    // ---- x value: read the leaf directly
    std::string leafName = "Oscillation_Parameters/#" +
                           std::to_string(lf->second.idx) +
                           "_" + par;
    auto* x_leaf = t->GetLeaf(leafName.c_str());
    if (!x_leaf) { std::cerr<<"Skip "<<stem<<" (leaf "<<leafName<<" not found)\n"; continue; }
    double xval = x_leaf->GetValue();

    std::cout<< stem.c_str() << ", " << leafName.c_str() << ", " << xval << ", " << llh << std::endl;

    pts[par].push_back({xval, llh});
  }

  if (pts.empty()) { std::cerr<<"No points collected\n"; return; }

  TFile fout("scan_curves_leaf.root","RECREATE");

  for (auto& [par,v] : pts)
  {
    std::sort(v.begin(), v.end(),
              [](auto& a, auto& b){ return a.first < b.first; });

    std::vector<double> x,y;
    for (auto& p : v){ x.push_back(p.first); y.push_back(p.second); }

    TGraph g(x.size(), x.data(), y.data());
    g.SetName((par+"_scan").c_str());
    fout.cd(); g.Write();

    const char* axis = leafInfo.at(par).label;            // ← lookup again
    TCanvas c;
    g.SetMarkerStyle(20); g.Draw("ALP");
    g.SetTitle((std::string(axis) + " scan;" +
               axis + "; totalLlhAtBestFit").c_str());        // ← use axis

    c.SaveAs((PLOT_DIR/(par+"_scan.pdf")).c_str());
  }


  fout.Close();
  std::cout<<"✓ plots → "<<PLOT_DIR<<"   graphs → scan_curves_leaf.root\n";
}
