# Oscillation Parameters Profiling
likelihood scan: freeze every other parameter at one fixed value  
profiling: fix one point and re-fit others

## 1) Make yaml files
```bash
g++ -std=c++17 -O2 make_osc_yaml_grid_profilling.cpp -lyaml-cpp -o make_osc_yaml_grid
./make_osc_yaml_grid  
```
## 2) Generate config files
Remember to disable oscillation parameters in the configParSet.yaml
```bash
sbatch config_grid_profiling.sh
```
## 3) Use GUNDAM
count JSON files
```bash
NUM=$(ls outputs/profiling_json_files/config_DUNE_With_*.json | wc -l)
echo "$NUM JSON files found"
```
submit one array task per file
```bash
sbatch --array=0-$((NUM-1)) config_grid_profiling_root.sh
```
## 4) Draw profiling results
```bash
root -l -b profiling_scan_plot.C
```
