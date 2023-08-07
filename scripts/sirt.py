import cv2 as cv
import mrcfile as mrc
import numpy as np
from tqdm import tqdm
import torch
def contrast(image):
    image = (image - np.min(image)) / (np.max(image) - np.min(image))
    t_sd = 40 / 255.0
    t_mean = 150/255.0
    mean = np.mean(image)
    sd = np.mean((image-mean)**2)
    sd = np.sqrt(sd)
    F = t_sd / sd
    A = t_mean - mean * t_sd / sd
    black = -A / F
    white = (1-A) / F
    image[image<black] = black
    image[image>white] = white
    image = image.astype(float)
    return image
# data_path = "/home/vrlab/limzh/dataset1/ihuman/Insulin_secretory_granule/20220105_data-05/20220105_data-05.rec/21"
# mrc_path = data_path + "/Position_21_1_ali.mrc"
# data_path = "/home/vrlab/limzh/dataset1/ihuman/neuron/20210929_neuron/20210929_rec/25"
# mrc_path = data_path + "/20210929_25_ali.mrc"
raw_path = "/home/vrlab/desktop/liufx/res-Etomo"
# data_path = "/home/vrlab/limzh/dataset1/ihuman/vesicles/20220121_grid5/20220121_grid5_rec/10"
# mrc_path = data_path + "/Position_1_2_ali.mrc"
mrc_path = raw_path + "/denoised_ali.mrc"
save_path = raw_path + "/denoised_ali_sirt.mrc"
images = mrc.open(mrc_path)
new_mrc = mrc.new(save_path,overwrite=True)
new_images = np.zeros(images.data.shape,dtype = np.float32)
print(images.data.shape)
print(np.max(images.data))
print(np.min(images.data))
images.print_header()
cutoff = 0.35
falloff = 0.035
for i in tqdm(range(images.data.shape[0])):
    # 1024 * 1024
    # image = np.log10((images.data[i]+1).astype(np.float32))
    image = images.data[i].copy()
    # image = cv.imread("lena.jpeg",cv.IMREAD_GRAYSCALE).astype(float)/255
    # print(image)
    # print(image.dtype)
    # fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(image),norm = 'ortho'))
    fft = torch.fft.rfft2(torch.from_numpy(image),norm='ortho').numpy()

    filter = np.zeros(fft.shape,dtype = np.float32)
    h,w = filter.shape
    cutoff = 4
    N = 30
    a = 0.00195
    maxfreq  = (2 ** 0.5) * max(h,w)
    vmax = int(maxfreq*cutoff)
    fmax = int(maxfreq*cutoff) / maxfreq
    for u in range(h):
        for v in range(w):
            if u < h//2:
                freq = np.sqrt(u**2+v**2)
                f = freq / maxfreq
            else:# [0,1,2,3,4,-5,-4,-3,-2,-1]
                freq = np.sqrt((u-h)**2+v**2)
                f = freq / maxfreq
            # f = v / w
            # filter[u,v] = np.exp(-f/(2 * cutoff ** 2))
            if freq < vmax:
                if f==0:
                    filter[u,v] = 0.2
                elif a / f >= 1:
                    filter[u,v] = 1
                elif a / f < 1:
                    filter[u,v] = (1-(1-a/f)**N) * freq
            else:
                
                # if u < h//2:
                #     old_u = int(u/((u*u+v*v)**0.5)*(int(maxfreq*cutoff)-1))
                #     old_v = int(v/((u*u+v*v)**0.5)*(int(maxfreq*cutoff)-1))
                #     freq = np.sqrt(old_u**2+old_v**2)
                #     arg = ((u*u+v*v)**0.5-int(maxfreq*cutoff))/max(0.01,maxfreq*falloff)
                # else:
                #     old_u = h-int((h-u)/(((h-u)*(h-u)+v*v)**0.5)*(int(maxfreq*cutoff)-1))-1
                #     old_v = int(v/(((h-u)*(h-u)+v*v)**0.5)*(int(maxfreq*cutoff)-1))
                #     freq = np.sqrt((h-old_u)**2+old_v**2)
                #     arg = (((h-u)*(h-u)+v*v)**0.5-int(maxfreq*cutoff))/max(0.01,maxfreq*falloff)
                
                # if filter[old_u,old_v] == 0:
                #     f = freq / maxfreq
                #     filter[old_u,old_v] = (1-(1-a/f)**N) * freq
                # filter[u,v] = filter[old_u,old_v] * (np.e **(-0.5*arg*arg))
                arg = (freq-vmax)/max(0.01,maxfreq*falloff)
                filter[u,v] = (1-(1-a/fmax)**N) * vmax * (np.e **(-0.5*arg*arg))
            
    filtered_fft = fft * filter
    filtered = torch.fft.irfft2(torch.from_numpy(filtered_fft),norm='ortho').numpy()
    # print(filtered.shape)
    new_images[i] = filtered
new_mrc.set_data(new_images)

# image = contrast(image)
# filtered = contrast(np.abs(filtered))
# cv.imshow('image',image)
# cv.waitKey()
# cv.imshow('filtered',np.abs(filtered))
# cv.waitKey()
# # cv.imshow('fft',np.log(1+np.abs(fft)))
# cv.imshow('fft',np.abs(fft))
# cv.waitKey()
# # cv.imshow('filtered_fft',np.log(1+np.abs(filtered_fft)))
# cv.imshow('filtered_fft',np.abs(filtered_fft))
# cv.waitKey()

