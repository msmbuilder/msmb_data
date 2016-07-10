# msmb_data

Testing and example data for MSMBuilder. The data is available under the
creative commons attribution license (CC-BY). Please cite:

 - fs_peptide http://dx.doi.org/10.6084/m9.figshare.1030363
 - met_enkephalin http://dx.doi.org/10.6084/m9.figshare.1026324
 - alanine_dipeptide http://dx.doi.org/10.6084/m9.figshare.1026131

# Installation

## `conda`

This package is designed to be used via conda

```bash
conda install msmb_data
```

## Manual installation

Copy the `msmb_data` directory to `/path/to/python/prefix/share/msmb_data`


# Using

You should probably use the methods in `msmbuilder` to load this data. To
do it by hand, try something like this

```python
import mdtraj
import sys
import glob
traj_fns = glob.glob("{}/share/msmb_data/alanine_dipeptide/trajectories-*.???"
                      .format(sys.prefix))
top_fn = "{}/share/msmb_data/alanine_dipeptide/alanine_dipeptide.pdb".format(sys.prefix)
trajs = [mdtraj.load(tfn, top=top_fn) for tfn in traj_fns]
```