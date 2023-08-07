import mrcfile
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

paths = glob.glob("/home/vrlab/desktop/liufx/mrc_files/rawdata_mrc/Position_10_*.mrc")
paths.sort(key=lambda x: float(x[60:63]))
# print(paths[0])

angle = []
for line in paths:
        angle.append(float(line.split("[")[1].split("]")[0]))
# print(angle)
angle = np.array(angle)
idx = angle.argsort()
sample = mrcfile.open(paths[0])

unaligned_mrc = np.zeros((len(paths),sample.data.shape[1],sample.data.shape[2]),dtype=np.float32)
for i,index in enumerate(idx):
    print(i)
    
    mrc_data = mrcfile.open(paths[index]).data
    # print(mrc_data.dtype)
    # mrc_data = np.sum(mrc_data.data.copy(),axis=0)
    # plt.imshow(mrc_data,cmap='gray')
    # plt.show()
    unaligned_mrc[i] = mrc_data
    # print(mrc_data.shape)

with mrcfile.new("/home/vrlab/desktop/liufx/mrc_files/rawdata_unaligned.mrc",overwrite=True) as mrc:
    mrc.set_data(unaligned_mrc)
    mrc.header.ispg = 0
    mrc.header.cella.x = 12404.055
    mrc.header.cella.y = 12404.055
    # mrc.header.cella.z = 11.178