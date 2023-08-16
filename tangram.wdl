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
        String docker = "izabellaleahz/tangramworkflow:latest"
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

        adata_sp = sc.read_h5ad("~{anndata_file_sp}")

        adata_sc = sc.read_h5ad("~{anndata_file_sc}")

        markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
        markers = list(np.unique(markers_df.melt().value.values))

        tg.pp_adatas(adata_sc, adata_sp, genes=markers)

        ad_map = tg.map_cells_to_space(adata_sc, adata_sp,
            mode="cells",
            density_prior='rna_count_based',
            num_epochs=500,
            device='cpu',
        )

        ad_map.write_h5ad('outputs/ad_map.h5ad')



        CODE
        gsutil -m rsync -r outputs ~{output_dir}
    >>>
    output {
        File tangram_object = 'outputs/ad_map.h5ad'
    }
    runtime { 
      docker: docker
      memory: memory
      bootDiskSizeGb: 12
      disks: "local-disk " + (ceil(size(anndata_file_sc, "GB")*4) + extra_disk_space) + " HDD"
      cpu: cpu
      preemptible: preemptible
      gpuType: "nvidia-tesla-t4-vws"
      gpuCount: 2
      zones: ["us-central1-c"] 
    }
}