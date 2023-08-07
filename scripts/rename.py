import os
import glob

ODD_paths = glob.glob("/media/vrlab/sdc1/denoise/training_data/partA/*.mrc")
for ODD_file in ODD_paths:
    new_ODD_file = ODD_file.replace('_ODD', '')
    os.rename(ODD_file, new_ODD_file)
print("finish renaming ODD files\n")

EVN_paths = glob.glob("/media/vrlab/sdc1/denoise/training_data/partB/*.mrc")
for EVN_file in EVN_paths:
    new_EVN_file = EVN_file.replace('_EVN', '')
    os.rename(EVN_file, new_EVN_file)
print("finish renaming EVN files\n")