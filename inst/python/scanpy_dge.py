import scanpy as sc

def scanpy_dge(adata, group_by, method = "wilcoxon"):
    """
    Runs differntial gene expression on an anndata object using Scanpy's
    rank_genes_groups function.
    
    Arguments
    ----------
    adata: Anndata object on which to run the test.
    
    groupby: the metadata variable used for comparing groups.
    
    method: statistical method to use for differential gene expression. The 
    default is "wilcoxon". 
    
    Returns
    ----------
    A pandas dataframe with DGE results. 
    """
    # Run rank_genes_groups on the object and the group by category
    sc.tl.rank_genes_groups(
        adata, 
        groupby = group_by, 
        method = method, 
        pts = True
        )
    
    # Fetch results in data frame format and return
    # Use group = None to return all entries
    results_df = sc.get.rank_genes_groups_df(r.adata, group = None)
    
    return results_df
