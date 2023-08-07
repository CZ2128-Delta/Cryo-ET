import numpy as np
import eer
import os
import mrcfile
import glob

eer_path = "/media/vrlab/sdc1/denoise/20230404_Ht22_Pos12_Pos17/Position_17/eer_data/*.eer"
data_path = "/media/vrlab/sdc1/denoise/20230404_Ht22_Pos12_Pos17/Position_17/"

paths = glob.glob(eer_path)
print(len(paths))
for path in paths:
    name = os.path.basename(path)[:-4] + ".mrc"

    pathA = data_path + "training_data/partA/" + name
    pathB = data_path + "training_data/partB/" + name
    pathC = data_path + "raw_data/" + name

    if os.path.exists(pathA):
        print(pathA)
        continue
    
    num_frames = eer.get_num_frames(eer.open_eer(path))
    mov = eer.read_eer(path,num_frames,1)
    mov_c = eer.read_eer(path,1,1)
    
    A = np.sum(mov[::2],axis=0).astype(np.float32)
    B = np.sum(mov[1:][::2],axis=0).astype(np.float32)
    
    with mrcfile.new(pathA,overwrite=True) as f:
        f.set_data(A)
    with mrcfile.new(pathB,overwrite=True) as f:
        f.set_data(B)
    with mrcfile.new(pathC,overwrite=True) as f:
        f.set_data(mov_c)