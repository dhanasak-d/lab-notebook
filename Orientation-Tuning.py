# %% [markdown]
# # Orientation Tuning


# %%
# general python modules for scientific analysis
import sys, pathlib, os
import numpy as np

sys.path += ['physion/src'] # add src code directory for physion
from physion.utils import plot_tools as pt

from physion.analysis.read_NWB import Data,\
    scan_folder_for_NWBfiles

dataset = scan_folder_for_NWBfiles(\
        os.path.join(os.path.expanduser('~'), 
            'DATA', 'Taddy', 'PN_shGrid1-2026'),
            #for_protocols=['tuning']
            )

# %%
dFoF_parameters = dict(\
    roi_to_neuropil_fluo_inclusion_factor=0., # no factor here
    neuropil_correction_factor = 0.7,
    method_for_F0 = 'sliding_percentile',
    percentile=5., # percent
    sliding_window = 5*60, # seconds
)


# %%
dataset
# %%
