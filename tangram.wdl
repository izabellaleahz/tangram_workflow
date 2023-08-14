version 1.0
workflow tangram {
    input {
        String output_directory
        File anndata_file_sc
        File anndata_file_sp
        
        #general parameters
        Int cpu = 24
        String memory = "128G"
        Int extra_disk_space = 32
        String docker = "izabellaleahz/tangram:latest"
        Int preemptible = 2
    }
    String output_directory_stripped = sub(output_directory, "/+$", "")
    call run_tangram {
        input:
            output_dir = output_directory_stripped,
            anndata_file_sc = anndata_file_sc,
            anndata_file_sp = anndata_file_sp,
            

            cpu=cpu,
            memory=memory,
            extra_disk_space = extra_disk_space,
            docker=docker,
            preemptible=preemptible
    }
    output {
        File tangram_object = run_tangram.tangram_object
    }
}
task run_tangram {
    input {
        String output_dir
        File anndata_file_sp
        File anndata_file_sc
        
        String memory
        Int extra_disk_space
        Int cpu
        String docker
        Int preemptible
    }
    command <<<
        set -e
        mkdir -p outputs
        python <<CODE

        import os, sys
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        import scanpy as sc
        import torch
        import tangram as tg

        adata_sp = sc.read_h5ad("../baysor_crc_CosMx.h5ad")

        adata_sc = sc.read_h5ad("../../scRNA/scRNA_CRC.h5ad")

        sc.pp.filter_genes(adata_sp, min_cells=1)
        sc.pp.filter_genes(adata_sc, min_cells=1)

        coords = np.array((adata_sp.obs["X"],adata_sp.obs["Y"]))
        adata_sp.obsm["spatial"] = coords.T 

        adata_sp.uns["spatial"] = adata_sp.obsm["spatial"]
        adata_sp.obs['total_counts'] = (adata_sp.X > 0).sum(axis=1)
        adata_sc.obs['total_counts'] = (adata_sc.X > 0).sum(axis=1)

        
        mapping_dict = {
            'B': 'B-cells/plasma cells',
            'T': 'T-cell',
            'Endothelial': 'endothelial (vascular)',
            'Epithelial': 'epithelial',
            'Fibroblast': 'fibroblast (ECM)',
            'Hepatocytes': 'hepatocyte',
            'Myeloid': 'macrophage/myeloid',
            'Mast': 'macrophage/myeloid',
            'Plasma': 'B-cells/plasma cells'
        }

        adata_sc.obs['cell_types_coarse'] = [mapping_dict[label] for label in adata_sc.obs['category']]

        mapping_dict = {
            2: 'P02',
            13: 'P13',
            18: 'P18',
            19: 'P19', 
            28: 'P28', 
            35: 'P35'
        }

        adata_sp.obs['patientID'] = [mapping_dict[label] for label in adata_sp.obs['sample']]

        sc.tl.rank_genes_groups(adata_sc, groupby="cell_types_coarse", use_raw=False)
        markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
        markers = list(np.unique(markers_df.melt().value.values))
        print("Number of markers: " + str(len(markers)))

        tg.pp_adatas(adata_sc, adata_sp, genes=markers)

        ad_map = tg.map_cells_to_space(adata_sc, adata_sp,
            mode="cells",
            # mode="clusters",
            # cluster_label='cell_types_coarse',  # .obs field w cell types
            density_prior='rna_count_based',
            num_epochs=500,
        #     device="cuda:0",
            device='cpu',
        )

        ad_map.write_h5ad('~{output_dir}' + 'ad_map.h5ad')



        CODE
        gsutil -m rsync -r outputs ~{output_dir}
    >>>
    output {
        File tangram_object = '~{output_dir}' + 'ad_map.h5ad'
    }
    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + (ceil(size(anndata_file_sc, "GB")*4) + extra_disk_space) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}