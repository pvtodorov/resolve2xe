# resolve2xe
convert Resolve biosciences outputs to xenium-like format that can be read by Seurat

to use:
1. `git clone https://github.com/pvtodorov/resolve2xe.git`
2. `cd resolve2xe`
3. `pip install -r requirements.txt` (preferably in a virtual environment of your choice)
4. call the script with
```
python resolve2xe.py --baysor_results_dir <path to baysor segmentation directory> \
                     --cellpose_roi_path <path to zipped Cellpose ROIs> \
                     --output_dir <path where you want to write Xenium-like output>
```
