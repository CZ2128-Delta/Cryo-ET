## ---------------------------------------------------------------------------
##    Copyright (c) 2021 Structura Biotechnology Inc. All rights reserved. 
##         Do not reproduce or redistribute, in whole or in part.
##      Use of this code is permitted only under licence from Structura.
##                   Contact us at info@structura.bio.
## ---------------------------------------------------------------------------

import numpy as n
import eerdecompressor

def open_eer(fname):
    # return eer object
    # this opens the file, reads it completely into RAM and reads all headers etc
    # by default, this will use b-splines and float output buffer
    return eerdecompressor.ElectronCountedFramesDecompressor(str(fname))
    
def read_eer_shape(fname):
    eer = open_eer(fname)
    nx, ny, nFrames = eer.getSize()
    return (nFrames, ny, nx)

def get_frame_shape(eer):
    # just get the frame shape
    nx, ny, nFrames = eer.getSize()
    return ny, nx # TODO CHECK THIS (okay when nx == ny)

def get_output_frame_shape(raw_frame_shape, upsampfactor):
    assert upsampfactor in [-32, -16, -8, -4, -2, 1, 2, 4]
    frame_shape = list((n.array(raw_frame_shape) * ( upsampfactor if upsampfactor > 0 else -1.0/upsampfactor )).astype(n.int))
    return frame_shape[0], frame_shape[1]

def get_num_frames(eer):
    # just get the num frames
    nx, ny, nFrames = eer.getSize()
    return nFrames

def read_eer_single_fraction(eer, startframe, endframe, upsampfactor=1, frame_buffer=None):
    # allocate buffer
    # fill buffer with data from requested frames
    # return buffer
    raw_frame_shape = get_frame_shape(eer)
    frame_shape = get_output_frame_shape(raw_frame_shape, upsampfactor)
    num_frames  = get_num_frames(eer)
    assert startframe < num_frames
    assert endframe <= num_frames
    if frame_buffer is None:
        # TODO this is uint8 for now but should be float32 for the Bspline version
        # frame_buffer = n.zeros(frame_shape, n.float32)
        frame_buffer = n.zeros(frame_shape, n.uint8)
    else:
        assert frame_buffer.shape == frame_shape
        # assert frame_buffer.dtype == n.float32
        assert frame_buffer.dtype == n.uint8
    for frameidx in range(startframe, endframe):
        # eer.decompressImage_AddTo_BSpline(frame_buffer, upsampfactor, frameidx)
        eer.decompressImage_AddTo(frame_buffer, upsampfactor, frameidx)
    return frame_buffer

def get_num_fractions(eer, num_frames_per_fraction, drop_trailing=True):
    num_frames = get_num_frames(eer)
    if drop_trailing:
        num_fractions = int(n.floor(float(num_frames) / float(num_frames_per_fraction))) # truncate / round down to drop trailing
    else:
        num_fractions = int(n.ceil(float(num_frames) / float(num_frames_per_fraction))) # round up to keep all
    return num_fractions

def get_num_frames_per_fraction(eer, num_fractions):
    num_frames = get_num_frames(eer)
    return int(n.floor(num_frames / float(num_fractions))) 

def read_eer_fractions(eer, num_frames_per_fraction, upsampfactor=1, drop_trailing=True):
    num_frames = get_num_frames(eer)
    num_fractions = get_num_fractions(eer, num_frames_per_fraction, drop_trailing)
    raw_frame_shape = get_frame_shape(eer)
    upsampfactor = int(upsampfactor)
    frame_shape = get_output_frame_shape(raw_frame_shape, upsampfactor)
    mov_shape = (num_fractions, ) + frame_shape
    mov = n.zeros(mov_shape, n.uint8)
    for frac_idx in range(num_fractions):
        startframe = frac_idx * num_frames_per_fraction
        endframe = min( (frac_idx + 1) * num_frames_per_fraction, num_frames)
        read_eer_single_fraction(eer, startframe, endframe, upsampfactor=upsampfactor, frame_buffer=mov[frac_idx])
    return mov

def read_eer(fname, num_fractions, upsampfactor):
    eer = open_eer(fname)
    num_frames_per_fraction = get_num_frames_per_fraction(eer, num_fractions)
    return read_eer_fractions(eer, num_frames_per_fraction, upsampfactor=upsampfactor, drop_trailing=True)

def read_eer_frames(fname, num_frames_required, upsampfactor):
    eer = open_eer(fname)
    num_frames = get_num_frames(eer)
    raw_frame_shape = get_frame_shape(eer)
    upsampfactor = int(upsampfactor)
    frame_shape = get_output_frame_shape(raw_frame_shape, upsampfactor)
    mov_shape = (num_frames_required, ) + frame_shape
    mov = n.zeros(mov_shape, n.uint8)
    idx_list = n.random.permutation(num_frames)
    for i in range(num_frames_required):
        frame_idx = idx_list[i]
        read_eer_single_fraction(eer, frame_idx, frame_idx+1, upsampfactor=upsampfactor, frame_buffer=mov[i])
    return mov

# https://github.com/fei-company/EerReaderLib
# https://github.com/cryoem-uoft/tfs_eer
# https://fei-box.app.box.com/s/3enc37s3nnfdh558ovysxus8bcyrc8t4/folder/117486307754

# Notes about how EER will work
# 
# at import time, set the number of frames and the upsamp factor
# at import time, assume user is providing psize as the raw psize of the 4k 
# but make sure to divide psize by upsamp factor
# 
# gain reference/defect file may come in 2 forms:
# TIFF file at 4096 x 4096, to be multiplied, with zeros for defects, must be re-scaled
# MRC file, same as above
# another source like empiar 10424 has to be converted to one of the above.
# 
# So with the gain:
# - at import time, make sure it's 4k x 4k
# - at use time, upscale to final upsamp size, splitting each value of 4, 16 etc positions 
# - at use time, find zeros in upscaled gain reference and those are defects (this can be generic not just EER)
# - if defect file also present, merge the defects
#
# Compilataion notes:
# - eerToolsNew folder from tfs_eer is all we need
# - makefile works (just run make clean && make)
# - eerdecompressor.so needs to be put in blobio dir
#
# Note: opening an EER file means reading the whole thing! only way to check length.
#