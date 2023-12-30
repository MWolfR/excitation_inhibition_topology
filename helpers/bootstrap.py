

def save_data(df, filename, data_id, dset):
    """
    Helper function to save results from the analysis of a single subvolume. 
    To be used to generate plots of mean +- std over subvolumes.
    """
    import os
    root = os.path.split(filename)[0]
    if not os.path.exists(root):
        os.makedirs(root)
    df.to_hdf(filename, key="{0}/{1}".format(dset, data_id))
