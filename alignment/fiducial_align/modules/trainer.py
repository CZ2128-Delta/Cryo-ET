import torch
import torch.nn as nn
import torch.functional as F

from typing import TYPE_CHECKING, NoReturn
if TYPE_CHECKING:
    from torch import Tensor
    from numpy import ndarray



class AlignTrainer:
    """ The trainer for fiducial alignment using torch ...
    """
    def __init__(self,
                 model,
                 optimizer, 
                 scheduler,
                 loss_fn,
                 fidpos_parser,
                 num_epochs: int = 40, 
                 **kwargs) -> None:
        self.model = model
        self.num_epochs = num_epochs
        self.loss_fn = loss_fn
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.fidpos_parser = fidpos_parser

        # TODO:
        


    def train_epoch(self, ) -> dict:
        # TODO: forward model and compute loss,  

        return {
            
        }
    
    def train(self, ) -> dict:
        # add scheduler
        for i in range(self.num_epochs):
            self.optimizer.zero_grad()
            uv_pred = self.model.forward_all_fidpos() # should be shape of (n, 2)
            loss = self.loss_fn(uv_pred, self.fidpos_parser.uv_gt)
            loss.backward()
            self.optimizer.step()

