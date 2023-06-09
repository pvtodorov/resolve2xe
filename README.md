# resolve2xe
convert Resolve biosciences outputs to xenium-like format that can be read by Seurat

to convert:
1. `git clone https://github.com/pvtodorov/resolve2xe.git`
2. `cd resolve2xe`
3. `pip install -r requirements.txt` (preferably in a virtual environment like [venv](https://docs.python.org/3/library/venv.html) or [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html))
4. call the script with
```bash
python resolve2xe.py --baysor_results_dir <path to baysor segmentation directory from Resolve> \
                     --cellpose_roi_path <path to zipped Cellpose ROIs from Resolve> \
                     --output_dir <path where you want to write Xenium-like output>
```

to load in R:
```r
source("LoadResolveBaysor.R")
resolve.obj = LoadResolveBaysor('<path to Xenium-like output>')
```
