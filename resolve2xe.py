#!/usr/bin/env python
# coding: utf-8


import argparse
import gzip
import requests
import pandas as pd
import sys
from pathlib import Path
import os
import csv
import glob
from collections import defaultdict
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
from shapely.geometry import MultiPoint, Polygon, Point
from rtree import index
import read_roi
from tqdm import tqdm

###  1.
###  to 10x matrix
###

def assign_indices(df):
    genes = defaultdict(lambda: len(genes))
    df['cell_idx'] = df['cell'].astype('category').cat.codes
    df['gene_idx'] = df['gene'].apply(lambda x: genes[x])
    return df, genes


def group_data(df):
    grouped_data = df.groupby(['cell_idx', 'gene_idx']).size().reset_index(name='counts')
    return grouped_data


def write_matrix(grouped_data, output_dir):
    rows, cols, data = grouped_data['cell_idx'].values, grouped_data['gene_idx'].values, grouped_data['counts'].values
    sparse_matrix = coo_matrix((data, (rows, cols)), shape=(max(rows) + 1, max(cols) + 1)).T
    with gzip.open(f"{output_dir}/cell_feature_matrix/matrix.mtx.gz", 'wb') as f:
        mmwrite(f, sparse_matrix, field='integer', precision=None, symmetry='general')


def query_api(genes):
    unique_genes = list(genes.keys())
    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
        json={
            'organism': 'mmusculus',
            'target': 'ENSG',
            'query': unique_genes,
        }
    )
    converted_genes = {gene: result['converted'] for gene, result in zip(unique_genes, r.json()['result'])}
    return converted_genes


def write_features(genes, converted_genes, output_dir):
    with gzip.open(f"{output_dir}/cell_feature_matrix/features.tsv.gz", 'wt') as f:
        for gene, idx in sorted(genes.items(), key=lambda x: x[1]):
            ensembl_id = converted_genes[gene]
            f.write(f"{ensembl_id}\t{gene}\tGene Expression\n")


def write_barcodes(df, output_dir):
    with gzip.open(f"{output_dir}/cell_feature_matrix/barcodes.tsv.gz", 'wt') as f:
        for cell in df['cell'].astype('category').cat.categories:
            f.write(f"{cell}\n")

#### 2.
#### takes in baysor_results.csv and generates convex hull boundaries for each cell
####


def make_roi_coordinates_df(roi_coordinates):
    print('Converting convex hulls to ROI table format.')
    coords_list = []
    for roi_name, coordinates in tqdm(roi_coordinates.items()):
        for x, y in coordinates:
            coords_list.append({'ROI_Name': roi_name, 'X': x, 'Y': y})
    coords_df = pd.DataFrame(coords_list)
    return coords_df


def baysor_results_df2coords_df(baysor_results_df):
    cell_groups = baysor_results_df.groupby('cell')
    convex_hulls = {}
    print('Calculating convex hulls for baysor segmentation.')
    for cell_name, group in tqdm(cell_groups):
        if cell_name != '0':  # Skip points not assigned to any cell
            points = group[['x', 'y']].values
            multi_point = MultiPoint(points)
            convex_hull = multi_point.convex_hull
            if convex_hull.geom_type == 'Polygon':
                coords = list(convex_hull.exterior.coords)
                convex_hulls[cell_name] = coords
    coords_df = make_roi_coordinates_df(convex_hulls)
    return coords_df

###
### 3
### CellPose roi_to_csv.py

def imagejzip2coords_df(zip_file_path):
    rois = read_roi.read_roi_zip(zip_file_path)
    rois_converted = []
    for roi_name, roi in rois.items():
        if roi['type'] == 'rectangle':
            x, y, w, h = roi['left'], roi['top'], roi['width'], roi['height']
            rois_converted.append({'ROI_Name': roi_name, 'X': x, 'Y': y})
        elif roi['type'] == 'freehand':
            x, y = roi['x'], roi['y']
            for i in range(len(x)):
                rois_converted.append({'ROI_Name': roi_name, 'X': x[i], 'Y': y[i]})
        else:
            raise ValueError(f'Unsupported ROI type: {roi["type"]}')
    coords_df = pd.DataFrame(rois_converted)
    return coords_df


def vertecies_to_shapely_polygon(roi_name=None, coordinates=None):
    if len(coordinates) < 3:
        print(f'Warning: Not enough coordinates for {roi_name} to form a valid polygon.')
        pg = Polygon()
    else:
        if coordinates[0] != coordinates[-1]:
            coordinates.append(coordinates[0])  # Add the first coordinate to the end if necessary
        pg = Polygon(coordinates)
    return pg


def make_polygon_df(coords_df):
    polygons = []
    for name, group in coords_df.groupby('ROI_Name'):
        vertices = [(x, y) for x, y in zip(group['X'], group['Y'])]
        # Add the last vertex identical to the first to make the polygon closed
        polygons.append({'ROI_Name': name, 'polygon': vertecies_to_shapely_polygon(name, vertices)})
    polygon_df = pd.DataFrame(polygons)
    polygon_df['Area'] = polygon_df['polygon'].apply(lambda x: x.area)
    polygon_df['valid'] = polygon_df['polygon'].apply(lambda x: x.is_valid)
    polygon_df['empty'] = polygon_df['polygon'].apply(lambda x: x.is_empty)
    polygon_df = polygon_df[polygon_df['valid']]
    polygon_df = polygon_df[~polygon_df['empty']]
    polygon_df = polygon_df.drop(['valid', 'empty'], axis=1)
    return polygon_df

###
###roi_containment.py
###


def map_nucleus2cell(nucleus_coords_df, cells_coords_df):
    print('Calculating nucleus to cell mappings.')
    nuc_cell_intersections = []

    idx = index.Index()
    cells_polygons = {}
    cells_names = []

    polygon_df_nuclei = make_polygon_df(nucleus_coords_df)
    polygon_df_cells = make_polygon_df(cells_coords_df)

    print('Building up rtree index of cells.')
    for cell_name, group in tqdm(polygon_df_cells.groupby('ROI_Name')):
#         print(cell_name, group)
        cell_polygon = group['polygon'].reset_index(drop=True)[0]
#         import pdb; pdb.set_trace()
        cells_polygons[cell_name] = cell_polygon
        cells_names.append(cell_name)
        idx.insert(len(cells_names) - 1, cell_polygon.bounds)

    print('Checking if nuclei are contained in cells.')
    for nucleus_name, group in tqdm(polygon_df_nuclei.groupby('ROI_Name')):
        nucleus_polygon = group['polygon'].reset_index(drop=True)[0]
        for i in idx.intersection(nucleus_polygon.bounds):
            cell_name = cells_names[i]
            cell_polygon = cells_polygons[cell_name]
            if cell_polygon.intersects(nucleus_polygon):
                nuc_cell_intersection = {'ROI_Name_nucleus': nucleus_name, 'ROI_Name_cell': cell_name}
                nuc_cell_intersections.append(nuc_cell_intersection)
    nuc_cell_df = pd.DataFrame(nuc_cell_intersections)
    return nuc_cell_df

###
###
### transcripts_in_nucleus.py


def get_transcripts_points(baysor_results_df):
    transcripts_points = [Point(x, y) for x, y in zip(baysor_results_df['x'], baysor_results_df['y'])]
    return transcripts_points


def create_rtree(polygon_df_nuclei):  
    """Create an R-tree spatial index for the nucleus boundaries."""
    idx = index.Index()
    for i, row in enumerate(polygon_df_nuclei.iterrows()):
        nucleus = row[1]['polygon']
        idx.insert(i, nucleus.bounds, obj=nucleus)
    return idx


def transcripts_in_nucleus(baysor_results_df, nucleus_coords_df):
    print('Loading nucleus boundaries...')
    nuclei = {}
    for nucleus_name, group in make_polygon_df(nucleus_coords_df).groupby('ROI_Name'):
        nuclei[nucleus_name] = group['polygon']
    print(f'Loaded {len(nuclei)} nuclei.')
    print('Loading transcripts...')
    transcripts = get_transcripts_points(baysor_results_df)
    print(f'Loaded {len(transcripts)} transcripts.')
    # Create R-tree spatial index for nuclei
    print('Creating R-tree index...')
    nucleus_polygon_df = make_polygon_df(nucleus_coords_df)
    idx = create_rtree(nucleus_polygon_df)
    # Check if each transcript is within a nucleus
    print('Checking transcripts...')
    results = []
    for i, transcript in tqdm(enumerate(transcripts), total=len(transcripts)):
        # Check if transcript is within any nucleus
        within_nucleus = False
        for jj in idx.intersection((transcript.x, transcript.y)):
            if transcript.within(nuclei[list(nuclei.keys())[jj]]).reset_index(drop=True)[0]:
                within_nucleus = True
                break
        results.append(within_nucleus)
    # Write results to output file
    baysor_results_df['within_nucleus'] = results
    return baysor_results_df



def make_cells_df(baysor_cell_stats, nucleus2cell_mapping, nucleus_coords_df):
    nucleus_polygon_df = make_polygon_df(nucleus_coords_df)
    # Apply transformations to cells DataFrame
    cells_df = (
        baysor_cell_stats.assign(
            control_probe_counts=0,
            control_codeword_counts=0,
            total_counts=lambda df: df['n_transcripts'],
            cell_area=lambda df: df['area'],
            x_centroid=lambda df: df['x'],
            y_centroid=lambda df: df['y'],
            cell_id=lambda df: df['cell']
        )
        .rename(columns={'n_transcripts': 'transcript_counts'})
    )
    cells_df = cells_df[['cell_id', 'x_centroid', 'y_centroid', 'transcript_counts', 'control_probe_counts', 'control_codeword_counts', 'total_counts', 'cell_area']]
    cells_df = pd.merge(cells_df, nucleus2cell_mapping, left_on='cell_id', right_on='ROI_Name_cell', how='left')
    cells_df = pd.merge(cells_df, nucleus_polygon_df, left_on='ROI_Name_nucleus', right_on='ROI_Name', how='left')
    cells_df = cells_df.rename(columns={"ROI_Name_nucleus": "nucleus_id", "Area": "nucleus_area"})
    cells_df = cells_df.drop(columns=['ROI_Name_cell', 'ROI_Name', "polygon"])
    return cells_df


def make_transcripts_df(baysor_results_df):
    transcripts_df = (
        baysor_results_df.assign(
            transcript_id = lambda df: df['molecule_id'],
            cell_id = lambda df: df['cell'],
            overlaps_nucleus = lambda df: df['within_nucleus'],
            feature_name = lambda df: df['gene'],
            x_location = lambda df: df['x'],
            y_location = lambda df: df['y'],
            z_location = 0,
            qv = 42.0
        )
    )
    transcripts_df = transcripts_df[['transcript_id', 'cell_id', 'overlaps_nucleus', 'feature_name', 'x_location', 'y_location', 'z_location', 'qv']]
    return transcripts_df


def make_cell_boundaries(coords_df):
    cell_boundaries_df = (
        coords_df.assign(
        cell_id = lambda df: df['ROI_Name'],
        vertex_x = lambda df: df['X'],
        vertex_y = lambda df: df['Y'],
        )
    )
    cell_boundaries_df = cell_boundaries_df[['cell_id', 'vertex_x', 'vertex_y']]
    return cell_boundaries_df


def convert_resolve_to_xenium(args):
    print('processing: ')
    for key, value in args.__dict__.items():
        print(key + ': ' + value)
    print('')
    # 1
    #create cell_feature_matrix outputs
    print("STAGE 1: Converting baysor results to 10x cell_feature_matrix.")  
    output_dir = args.output_dir
    cell_feature_matrix_dirpath = Path(f"{output_dir}/cell_feature_matrix/")
    cell_feature_matrix_dirpath.mkdir(exist_ok=True, parents=True)
    baysor_results_dir = args.baysor_results_dir
    baysor_results_path = f"{baysor_results_dir}/baysor_results.csv"
    baysor_results_df = pd.read_csv(baysor_results_path)
    baysor_results_df, genes = assign_indices(baysor_results_df)
    grouped_data = group_data(baysor_results_df)
    write_matrix(grouped_data, output_dir)
    converted_genes = query_api(genes)
    write_features(genes, converted_genes, args.output_dir)
    write_barcodes(baysor_results_df, args.output_dir)
    print(f"Matrix dimensions: {grouped_data.shape}")
    print(f"Number of genes: {len(genes)}")
    print(f"Number of cells: {len(baysor_results_df['cell'].astype('category').cat.categories)}")
    print("done.")
    print("")

    # 2
    # generate convex hulls for baysor segmentations
    print('STAGE 2: Generating cell segmentations from Baysor transcript assignment.')
    cells_coords_df = baysor_results_df2coords_df(baysor_results_df)
    cell_boundaries_df = make_cell_boundaries(cells_coords_df)
    print('writing cell_boundaries.csv.gz')
    cell_boundaries_df.to_csv(f"{output_dir}/cell_boundaries.csv.gz", index=False, compression='gzip')
    print("done.")
    print("")
    
    # 3
    #
    print('STAGE 3: Converting ImageJ nucleus ROIs.')
    cellpose_roi_path = args.cellpose_roi_path
    nucleus_coords_df = imagejzip2coords_df(cellpose_roi_path)
    print("done.")
    print("")
    
    # 4
    #
    print('STAGE 4: Mapping nuclei to cells.')
    nucleus2cell_mapping = map_nucleus2cell(nucleus_coords_df, cells_coords_df)
    print("done.")
    print("")
    
    # 5
    #
    print("STAGE 5: Checking for transcripts in nuclei.")
    baysor_results_df = transcripts_in_nucleus(baysor_results_df, nucleus_coords_df)
    transcripts_df = make_transcripts_df(baysor_results_df)
    print("writing transcripts.csv.gz")
    transcripts_df.to_csv(f"{output_dir}/transcripts.csv.gz", index=False, compression='gzip')
    print("done.")
    print("")
    
    # 6
    #
    print("STAGE 6: Writing cells.csv.gz")
    baysor_cell_stats_path = f"{baysor_results_dir}/baysor_cell_stats.csv"
    baysor_cell_stats = pd.read_csv(baysor_cell_stats_path)
    cells_df = make_cells_df(baysor_cell_stats, nucleus2cell_mapping, nucleus_coords_df)
    cells_df.to_csv(f"{output_dir}/cells.csv.gz", index=False, compression='gzip')
    print("done.")
    print("")


####
###
##
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process input file to create a cell feature matrix.')
    parser.add_argument('baysor_results_dir', type=str, help='path to baysor results directory')
    parser.add_argument('cellpose_roi_path', type=str, help='path to cellpose ROI zip')
    parser.add_argument('output_dir', type=str, help='path to output folder')
    args = parser.parse_args()
    convert_resolve_to_xenium(args)


