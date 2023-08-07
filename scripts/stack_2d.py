import mrcfile
import numpy as np
import glob

paths = glob.glob("/media/vrlab/sdc1/denoise/aligned_raw_data/denoised_results/20210929_26_ali_*.mrc")
paths.sort(key=lambda x: float(x[76:78]))

# paths_have_angle = glob.glob("/media/vrlab/sdc1/20210929_neuron/20210929_frames/20210929_1_*.eer")
# paths_have_angle.sort(key=lambda x: float(x[61:64]))
# print(paths)
# print(paths_have_angle)

#------ motioncor has sorted views by angle ------#

# angle = []
# for line in paths_have_angle:
#         angle.append(float(line.split("[")[1].split("]")[0]))
# # print(angle)
# angle = np.array(angle)
# idx = angle.argsort()
# # print(idx)
sample = mrcfile.open(paths[0])

unaligned_mrc = np.zeros((len(paths),sample.data.shape[0],sample.data.shape[1]),dtype=np.float32)
for i,index in enumerate(paths):
    print(i)
    
    mrc_data = mrcfile.open(paths[i]).data
    # print(mrc_data.dtype)
    # mrc_data = np.sum(mrc_data.data.copy(),axis=0)
    # plt.imshow(mrc_data,cmap='gray')
    # plt.show()
    unaligned_mrc[i] = mrc_data
    # print(mrc_data.shape)

with mrcfile.new("/media/vrlab/sdc1/reconstruct/Etomo/20210929_26_ali_denoised/20210929_26_ali_denoised.mrc",overwrite=True) as mrc:
    mrc.set_data(unaligned_mrc)
    mrc.header.ispg = 0
    mrc.header.cella.x = 12404.055
    mrc.header.cella.y = 12404.055
    # mrc.header.cella.z = 11.178