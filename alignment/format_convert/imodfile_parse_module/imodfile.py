import os
import sys
import json
import pickle
import numpy as np
from pathlib import Path
from imodfile_parse_module.constants import const 
from imodfile_parse_module.utils import *

from chunk import Chunk # a small toolkit for parsing IFF format

from typing import Union, List, Dict, TYPE_CHECKING
if TYPE_CHECKING:
    ...



# IMOD Binary File Format V1.2
# IMOD binary file format is similar to IFF format standard that uses chunk id's for data headings. Each chunk is 4 bytes long and is defined as a string of 4 chars. All numbers are stored in big-endian format regardless of the machine architecture.

class IMODFile(object):
    def __init__(self, file_path: Union[Path, str]) -> None:
        # The header's field '' defines the # of objts a `.mod` file contains
        self.header: Dict = None
        self.objts: List[Dict] = [] # each item is a dict of an `object`
        self.conts: List[Dict] = []
        self.minx: Dict = None

        self.file_path: Path = Path(file_path)
        self.suffix: str = self.file_path.suffix
        
        assert self.suffix in const.IMOD_SUPPORTED_FORMATS, f" *{self.suffix} format is not supported in IMOD support list: {const.IMOD_SUPPORTED_FORMATS}"
        # TODO: Here we only handle Vector Format files whose suffcies are .mod / .fid
        # Such file is the input format of `tiltalign` program.

        fstream = self.file_path.open('rb')
        # TODO: To add sanity-check
        while True:
            chunk_id = fstream.read(4) # read 4 bytes ID
            # determine which to parse ...
            if chunk_id == const.IMOD_ID:
                self.header = parse_header(fstream)
            elif chunk_id == const.OBJT_ID:
                self.objts.append(parse_objt(fstream))
            elif chunk_id == const.CONT_ID:
                self.conts.append(parse_cnts(fstream))
            elif chunk_id == const.MINX_ID:
                self.minx = parse_minx(fstream)
            elif chunk_id == const.IEOF_ID:
                break
            else: # undefined ID, we just skip the extra chunk.
                chunk_size = np.frombuffer(fstream.read(4), np.dtype('>i4'))[0]
                fstream.read(chunk_size)
    
    def writeContourData(self, track_path: str, bin: int) -> None:
        # Read tracking fiducial json file
        with open(track_path) as f:
            track_data = json.load(f)
        num_conts = 0
        contour_data = {}
        for i in range(track_data['frame_num']):
            num_conts = 4 # max(num_conts, len(track_data[f'frame_{i}'].keys()))
            cont_keys = track_data[f'frame_{i}'].keys()
            # # remove bad-tracked fiducials
            # cont_keys = list(cont_keys)
            # abondoned_conts = [3, 4, 8, 9, 18, 19, 20, 22, 23, 26, 27, 28, 30, 31, 33, 34, 35, 36, 37, 38]
            # for j in abondoned_conts:
            #     if f'group_{j}' in cont_keys:
            #         cont_keys.remove(f'group_{j}')
            # print(cont_keys)
            # ##############################
            # choose 4 best fiducials
            cont_keys = ['group_6', 'group_13', 'group_17', 'group_21']
            for key in cont_keys:
                if key not in contour_data:
                    contour_data[key] = []
                track_data[f'frame_{i}'][key][0] *= bin
                track_data[f'frame_{i}'][key][1] *= bin
                track_data[f'frame_{i}'][key].append(float(i))
                contour_data[key].append(track_data[f'frame_{i}'][key]) # [x, y, z]
        
        # Write tracking fiducial into new IMOD binary file
        new_data = []
        fstream = self.file_path.open('rb')
        while True:
            chunk_id = fstream.read(4) # read 4 bytes ID
            new_data.append(chunk_id)
            # determine which to parse ...
            if chunk_id == const.IMOD_ID: # write header chunk
                new_data.append(write_header(fstream))
            elif chunk_id == const.OBJT_ID:
                # write new object data
                # TODO: now num(object) = 1, need to consider num(object) > 1
                orig_cont_size, objt_info = write_objt(fstream, num_conts)
                new_data.extend(objt_info)
            elif chunk_id == const.CONT_ID:
                # skip original contour data
                flags_time_surf = skip_cnts(fstream)
                for i in range(1, orig_cont_size):
                    cont_chunk_id = fstream.read(4)
                    skip_cnts(fstream)
                # write new contour data
                first_cont = True
                for cont_num in [6, 13, 17, 21]: # range(1, num_conts + 1):
                    # # remove bad-tracked fiducials
                    # if cont_num in abondoned_conts:
                    #     continue
                    # ##############################
                    if first_cont == False:
                        new_data.append(cont_chunk_id)
                    first_cont = False
                    new_data.append(write_cnts(cont_num, contour_data, flags_time_surf))
            elif chunk_id == const.IEOF_ID:
                break
            else: # write extra chunk
                chunk_size = np.frombuffer(fstream.read(4), np.dtype('>i4'))[0]
                binary_chunk_size = struct.pack('>i', chunk_size)
                new_data.append(binary_chunk_size)
                chunk = fstream.read(chunk_size)
                new_data.append(chunk)
        binary_new_data = b''.join(new_data)
        print(binary_new_data)
        save_data_path = Path(f'{self.file_path.parent}/{self.file_path.stem}_track.fid')
        data_writer = save_data_path.open('wb')
        data_writer.write(binary_new_data)
        data_writer.close()
