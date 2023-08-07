#!/bin/bash

topaz denoise -m /media/vrlab/sdc1/denoise/20230404_Ht22_Pos12_Pos17/Position_17/saved_models/model_epoch100.sav \
              --patch-size 1024 \
              -o /media/vrlab/sdc1/denoise/HIV/denoised_results \
              /media/vrlab/sdc1/denoise/HIV/undenoised_views/TS_01_*.mrc