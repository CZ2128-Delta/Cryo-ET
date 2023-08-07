import numpy as np

import torch
import torch.nn as nn
import torch.functional as F
import os, sys

import mrcfile
from utils import *
from modules import AlignModel

import argparse

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from numpy import ndarray
    from torch import Tensor

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # training parameters
    parser.add_argument('--device', type=str, default='cuda', choices=['cpu', 'cuda'], help=' Device on which we optimize unkowns ...')
    parser.add_argument('--lr', type=float, default=1e-3, help=" Training learning rate for optimizer ...")
    # file path variables
    parser.add_argument('--prealign_path', type=str, default='/home/Assets/tilt_data/01_preali.mrc')
    parser.add_argument('--fidpos_path', type=str, default='/home/Assets/tilt_data/fidpos')
    parser.add_argument('--aligned_path', type=str, default='/home/Assets/tilt_data/01_ali.mrc') # dev only
    # entries to program tiltalign
    # basic info
    parser.add_argument('--unbinned_pix_size', type=float, default=.0675, help=" Nanometer/pixel, representing spatial resolution ...")
    parser.add_argument('--images_binned_num', type=int, default=1, help=' Bin number of the images ...')
    
    # ########## Tilt angle and rotation angle related options ###########
    # uv rot along z-axis
    parser.add_argument('--rot_option', type=int, default=1, choices=[0, 1, 2, 3, 4, -1], help=(" Type of rotation solution." 
                                                                                                "(0: for all rotations fixed at the initial angle;"
                                                                                                "1: for each view having an indepent rotation;"
                                                                                                "2: to enter general mapping of rotation variables"
                                                                                                "3 or 4: fir automapping of rotation variables"
                                                                                                "-1: to solve for a single rotation variable."))

    parser.add_argument('--init_rot', type=float, default=84.1, help=" Initial rotation angle along z-axis, to be optimized ...")
    # tilt related options
    parser.add_argument('--y_tilt_first', type=float, default=-57.0, help=" The tilt angle of the first view, in degrees. Use this option together with 'tilt_incr' ...")
    parser.add_argument('--y_tilt_incr', type=float, default=3.0, help=" Increment between tilt angles ...")
    parser.add_argument('--y_tilt_option', type=int, default=1, help=(" Type of tilt angle solution:"
                "0 to fix all tilt angles at their initial values,"
                "1 to solve for all tilt angles except for a specified view,"
                "2 to solve for all tilt angles except for the view at minimum tilt,"
                "3 to solve for all tilt angles except for a specified view and the view at minimum tilt,"
                "4 to specify a mapping of tilt angle variables,"
                "5 or 6 to automap groups of tilt angles (5 for linearly changing values or 6 for values all the same within a group), or"
                "7 or 8 to automap and fix two tilt angles (7 for linearly cchanging values or 8 for values all the same within a group)")
                )

    parser.add_argument('--is_robust_fitting', type=bool, default=True, help=" Whether use robust fitting refining results or not ...")
    parser.add_argument('--init_xyz_path', type=str, default='...', help=" Path of initial xyz of fidicuial points ...")
    # TODO: Please add more options input if needed
    # ...

    args = parser.parse_args()

    # TODO: parse input fid_pos file from imod 3d-mod file ... 
    fidpos_parser = FidPosParser(file_path=args.fidpos_path) 
    uv_list = fidpos_parser.get_all_uv()
    uv_list = torch.from_numpy(uv_list).cuda()

    # TODO: create alignment model 
    model = AlignModel().to(args.device) 
    uv_pred_list = model.forward_all_fidpos()
    diff = uv_list - uv_pred_list
    # print(uv_pred_list[0], uv_list[0])

    print(model.forward_one_fidpos(14, 0, True))

    # TODO: create optimizer for align model
    optimizer = torch.optim.Adam([
        {'params':model.parameters(), 'lr':args.lr}
    ], lr=args.lr
    )

    # TODO: trainer ...
    fid_align_trainer = ...


