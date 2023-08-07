import numpy as np
import mrcfile

from typing import Union, Any, TYPE_CHECKING, NoReturn, Optional
if TYPE_CHECKING:
    from numpy import ndarray

# Centroid positions of ficuial markers are stored in .pkmod file, which can be converted into .wimp format.
# This is a customed parser for .wimp to read imod uv positions out ... 

class FidPosParser(object):
    """ Parser for fiducial marker positions.
    """
    def __init__(self, file_path: Optional[str]=None) -> None:
        self.is_file_loaded = self.load_file(file_path)
    
    def load_file(self, file_path: Optional[str]=None) -> bool:
        """ Load and parse wimp file. 
        """
        if file_path is None:
            return False

        # TODO: ########### Parse file ##########
        # 
        # 
        self.markers_uv = dict() # we store marker positions using a dict -- {idx of marker: {ith view: position}}
        file = open(file_path, 'r')
        reading_view = False
        for line in file.readlines():
            if line == "\n": # reach the end of the file
                reading_view = False
                break
            if line[0:10] == "  Object #": # read a new marker
                reading_view = False
                j = line.split(':')[1].strip()
                self.markers_uv[f"marker_{j}"] = dict()
            if line[0:27] == "     #    X       Y       Z": # read the views of current marker
                reading_view = True # the following lines are views
                continue
            if reading_view:
                data = line.split()
                u, v, i = float(data[1]), float(data[2]), int(float(data[3]))
                self.markers_uv[f"marker_{j}"][f"view_{i}"] = np.array((u,v))
        #
        # #######################################
        return True
            
    def detach_file(self,) -> NoReturn:
        """ Detach loaded file.
        """
        self.is_file_loaded = False
    
    def normalize_uv(self, ) -> dict:
        # calculate u_max, v_max, u_min, v_min
        num_markers = len(self.markers_uv)
        views_list = list(self.markers_uv.values())
        uv_list = list()
        for j in range(num_markers):
            uv_list += list(views_list[j].values())
        uv_ndarray = np.array(uv_list)
        uv_max = np.amax(uv_ndarray, axis=0)
        uv_min = np.amin(uv_ndarray, axis=0)

        # normalize u,v
        normalized_uv = self.markers_uv
        for marker_key in normalized_uv.keys():
            for view_key in normalized_uv[marker_key].keys():
                normalized_uv[marker_key][view_key] -= uv_min
                normalized_uv[marker_key][view_key] /= uv_max - uv_min
        return normalized_uv

    def get_all_uv(self, ) -> 'ndarray':
        num_markers = len(self.markers_uv)
        views_list = list(self.markers_uv.values())
        uv_list = list()
        for j in range(num_markers):
            uv_list += list(views_list[j].values())
        uv_ndarray = np.array(uv_list)
        return uv_ndarray

    def get_uv(self, i: int, j: int) -> 'ndarray':
        """ Get j_th marker position of i_th view ... 

        Args:
            i (int): index of view 
            j (int): index of marker

        Returns:
            ndarray: uv positions
        """

        if not self.is_file_loaded:
            print("Warning: no file loaded, return (None, None)")
            return (None, None)
        
        marker_key = f"marker_{j}"
        if marker_key not in self.markers_uv.keys():
            # TODO: if not, raise a warning and return (None, None)
            ...
        
        view_key = f"view_{i}"
        if view_key not in self.markers_uv[marker_key].keys():
            # TODO: if not, raise a warning and return (None, None)
            ...

        return self.markers_uv[view_key][marker_key]

