import torch
import torch.nn as nn
import torch.functional as F

import numpy as np

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from numpy import ndarray
    from torch import Tensor


class AlignModel(nn.Module):
    """ A model containing alignment unknowns to be optimized via Adam. 
    """

    def __init__(self, **kwargs) -> None:
        super().__init__()
        # optimizable variables should store in nn.Parameter ... 
        # before this, should calculate accurate number of variables based on mapping and grouping strategy ... 
        #  
        # FIXME: Not for now, should fix it when HIV case works.
        # (1) here not considering distortion (skew and so on)
        # (2) hardcoded ... 

        file1 = open("/media/vrlab/sdc1/alignment/align.log", 'r')
        reading = False
        rotation_list = []
        mag_list = []
        tilt_list = []
        for line in file1.readlines():
            if line[0:16] == " view   rotation":
                reading = True
                continue
            if line == ' \n' and reading == True:
                reading = False
                break
            if reading == True:
                data = line.split()
                rotation_list.append(float(data[1]))
                mag_list.append(float(data[4]))
                if int(data[0]) in [1, 4, 9, 21, 27, 33, 37, 40]:
                    tilt_list.append(float(data[2]))
        rotation = np.array(rotation_list)
        mag = np.array(mag_list)
        tilt = np.array(tilt_list)

        file2 = open("/media/vrlab/sdc1/alignment/01.tltxf")
        offset_list = []
        for line in file2.readlines():
            data = line.split()
            offset_list.append(float(data[4]))
            offset_list.append(float(data[5]))
        offset = np.array(offset_list).reshape((40,2))

        self.ali_variables = nn.ParameterDict({
            'rot_angles': nn.Parameter(torch.from_numpy(rotation).float()), # FIXME hardcoding for now
            'mag': nn.Parameter(torch.from_numpy(mag).float()),  
            'xyz': nn.Parameter(torch.Tensor([[3025.53,1901.81,204.56], [2267.39,2223.36,30.38], [1832.58,1392.87,500.66], [1624.62,2982.10,793.48]])), #FIXME: load from input please. 3D coordinates of fiducial markers, # of vars: 3 * num of fid_markers
            # 'xyz': nn.Parameter(torch.Tensor([[ 838.00, -223.23, -177.71], [79.86, 98.32, -351.89], [-354.95, -732.16, 118.39], [-562.91, 857.07, 411.21]])), #FIXME: load from input please. 3D coordinates of fiducial markers, # of vars: 3 * num of fid_markers
          
            # 'tilt_angles': nn.Parameter(torch.Tensor([-57 + 3*(i-1) for i in [1, 4, 9, 21, 27, 33, 37, 40]])), # 
            'tilt_angles': nn.Parameter(torch.from_numpy(tilt).float()), # 
            'offset': nn.Parameter(torch.from_numpy(offset).float())
        }) 
        self.n_views = 40 # hardcoded
        self.n_markers = 4 # hardcoded
        self.device = 'cuda'
    

    def forward_one_fidpos(self, i: int, j: int, return_all: bool=False) -> 'Tensor':

        # TODO: sanity-check: if jth marker does not exist in ith view, we just return warning and return None
        ...
        

        # generate R_i
        rot_angle = self.ali_variables['rot_angles'][i] / 180 * torch.pi
        rot_matrix = torch.vstack([
            torch.hstack((torch.cos(rot_angle), -torch.sin(rot_angle))), 
            torch.hstack((torch.sin(rot_angle), torch.cos(rot_angle)))
            ],) 

        # generate Y_i
        check_list = np.array([0, 3, 8, 20, 26, 32, 36, 39])
        idx_1 = np.where(check_list <= i)[0][-1]
        idx_2 = idx_1 + 1
        
        if i == check_list[idx_1]:
            tilt_angle = self.ali_variables['tilt_angles'][idx_1] / 180 * torch.pi
        elif i == 14:
            tilt_angle = torch.tensor(-15.0, device=self.device) / 180 * torch.pi
        else:

            angle_idx_1 = self.ali_variables['tilt_angles'][idx_1].detach()  / 180 * torch.pi
            angle_idx_2 = self.ali_variables['tilt_angles'][idx_2].detach()  / 180 * torch.pi

            tilt_angle = angle_idx_1 + (angle_idx_2 - angle_idx_1) / (check_list[idx_2] - check_list[idx_1]) * (i-check_list[idx_1]) 

        tilt_matrix = torch.vstack([
            torch.hstack([torch.cos(tilt_angle), torch.tensor([0.0], device=self.device), torch.sin(tilt_angle)]),
            torch.tensor([0.0, 1.0, 0.0], device=self.device),
        ])

        # generate D_i 
        if i != 0:
            mag = self.ali_variables['mag'][i]
        else:
            mag = torch.tensor([1.0], device=self.device)
        
        distort_matrix = torch.vstack([
            torch.hstack([mag, torch.tensor([0.0, 0.0], device=self.device)]),
            torch.hstack([torch.tensor([0.0], device=self.device), mag, torch.tensor([0.0], device=self.device)]),
            torch.hstack([torch.tensor([0.0, 0.0], device=self.device), mag]),
        ])

        # generate xyz
        xyz = self.ali_variables['xyz'][j] # shape of (3, )
        xyz_hat = torch.mean(self.ali_variables['xyz'], dim=0)
        print(xyz_hat)
        # generate displacement vector
        if i != 0:
            offset = self.ali_variables['offset'][i]
        else:
            offset = torch.tensor([0.0, 0.0], device=self.device)
        A = (rot_matrix @ (tilt_matrix)) @ (distort_matrix)
        
        uv_pred = A @ (xyz) + offset 

        if return_all:
            return {
                'rot_angle': rot_angle.data * 180 / np.pi,
                'xyz': xyz.data,
                'mag': mag.data,
                'tilt_angle': tilt_angle.data * 180 / np.pi,
                'offset': offset.data,
                'uv': uv_pred.data,
                '222': rot_matrix
            }
        return uv_pred

    def forward_all_fidpos(self, ) -> 'Tensor':
        uv_pred_list = []
        for j in range(self.n_markers):
            for i in range(self.n_views):
                if j == 1:
                    if i < 18:
                        continue
                uv_pred_list.append(self.forward_one_fidpos(i, j))
        return torch.vstack(uv_pred_list)

if __name__ == '__main__':
    model = AlignModel().cuda()
