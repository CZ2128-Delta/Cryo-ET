#!/bin/bash

topaz denoise -a /media/vrlab/sdc1/denoise/20230404_Ht22_Pos12_Pos17/Position_17/training_data/partA \
              -b /media/vrlab/sdc1/denoise/20230404_Ht22_Pos12_Pos17/Position_17/training_data/partB \
              --save-prefix /media/vrlab/sdc1/denoise/20230404_Ht22_Pos12_Pos17/Position_17/saved_models/model \
              --num-workers 2