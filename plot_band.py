"""
Below is a quick script to show how to customize the plot
More functions can be found at https://phonopy.github.io/phonopy/phonopy-module.html
"""

import phonopy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

phonon = phonopy.load("example-Si8/phonopy_disp.yaml",
                      force_constants_filename = "example-Si8/FORCE_CONSTANTS")
phonon.auto_band_structure()
n = len([x for x in phonon._band_structure.path_connections if not x])
fig = plt.figure()
axs = ImageGrid(fig,
                111,  # similar to subplot(111)
                nrows_ncols=(1, n),
                axes_pad=0.11,
                label_mode="L")
phonon._band_structure.plot(axs)
for ax in axs:
    for line in ax.get_lines():
        line.set_color('black')
        #ax.set_ylim([0, 35])
plt.savefig('Test_band1.png')
