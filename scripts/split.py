import mrcfile
import numpy as np

stack_path = "/media/vrlab/sdc1/denoise/HIV/01.mrc"
views_path = "/media/vrlab/sdc1/denoise/HIV/undenoised_views/TS_01_"

stack = mrcfile.open(stack_path)
mic = np.array(stack.data, copy=False).astype(np.float32)
for i in range(mic.shape[0]):
    print(views_path + str(i+1) + ".mrc")
    with mrcfile.new(views_path + str(i+1) + ".mrc", overwrite=True) as view:
        view.set_data(mic[i])
        view.header.ispg = 0
        view.header.cella.x = 5008.5
        view.header.cella.y = 5181.3