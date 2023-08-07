from pathlib import Path
from imodfile_parse_module.constants import const
from typing import TYPE_CHECKING, Dict
import numpy as np
import struct
if TYPE_CHECKING:
    from io import BytesIO

def parse_header(fstream:'BytesIO') -> Dict:
    """ Parse `model` structure data. When encountering with `IMOD` ID, this function will be invoked to parse the following chunk data.
        There are totally 232 bytes in a header of standard `.mod` file.
    """
    version = fstream.read(4).decode('utf-8')
    name = fstream.read(128).decode('utf-8')
    x_max, y_max, z_max = np.frombuffer(fstream.read(12), np.dtype('>i4'))
    obj_size = np.frombuffer(fstream.read(4), np.dtype('>i4'))
    flags = np.frombuffer(fstream.read(4), np.dtype('>u4'))
    draw_mode, mouse_mode, black_lvl, white_lvl = np.frombuffer(fstream.read(16), np.dtype('>i4'))
    x_offset, y_offset, z_offset = np.frombuffer(fstream.read(12), np.dtype('>f4'))
    x_scale, y_scale, z_scale = np.frombuffer(fstream.read(12), np.dtype('>f4'))
    obj, cnt, pnt = np.frombuffer(fstream.read(12), np.dtype('>i4'))
    res, thresh = np.frombuffer(fstream.read(8), np.dtype('>i4'))
    pix_size = np.frombuffer(fstream.read(4), np.dtype('>f4'))
    units, csum = np.frombuffer(fstream.read(8), np.dtype('>i4'))
    alpha, beta, gamma = np.frombuffer(fstream.read(12), np.dtype('>f4'))

    return {
        'version': version,
        'name': name,
        'x_max': x_max,
        'y_max': y_max,
        'z_max': z_max,
        'obj_size': obj_size,
        'flags': flags,
        'draw_mode': draw_mode,
        'mouse_mode': mouse_mode,
        'blk_lvl': black_lvl,
        'wht_lvl': white_lvl,
        'x_offset': x_offset,
        'y_offset': y_offset,
        'z_offset': z_offset,
        'x_scale': x_scale,
        'y_scale': y_scale,
        'z_scale': z_scale,
        'obj_cur': obj,
        'cnt_cur': cnt,
        'pnt_cur': pnt,
        'res': res,
        'thresh': thresh,
        'pix_size': pix_size,
        'units': units,
        'csum': csum,
        'alpha': alpha,
        'beta': beta,
        'gamma': gamma,
    }

def parse_objt(fstream: 'BytesIO') -> Dict:
    name = np.frombuffer(fstream.read(64), np.byte)
    extra_data = np.frombuffer(fstream.read(64), np.dtype('>u4'))
    cont_size = np.frombuffer(fstream.read(4), np.dtype('>i4'))[0]
    flags = np.frombuffer(fstream.read(4), np.dtype('>u4'))[0]
    axis, draw_mode = np.frombuffer(fstream.read(8), np.dtype('>i4'))
    red, green, blue = np.frombuffer(fstream.read(12), np.dtype('>f4'))
    pdraw_size = np.frombuffer(fstream.read(4), np.dtype('>i4'))[0]
    symbol, symsize, lw2D, lw3D, linesty, symflags, sympad, trans = np.frombuffer(fstream.read(8), np.ubyte)
    mesh_size, surf_size = np.frombuffer(fstream.read(8), np.dtype('>i4'))

    return {
        'name': name,
        'extra': extra_data,
        'cont_size': cont_size, 
        'flags': flags, 
        'axis': axis,
        'draw_mode': draw_mode, 
        'red': red, 
        'green': green, 
        'blue': blue, 
        'pdraw_size': pdraw_size, 
        'symbol': symbol, 
        'sym_size': symsize,
        'line_width_2D': lw2D, 
        'line_width_3D': lw3D, 
        'linesty': linesty, 
        'symflags': symflags, 
        'sympad': sympad, 
        'trans': trans, 
        'mesh_size': mesh_size, 
        'surf_size': surf_size, 
    }

def parse_cnts(fstream: 'BytesIO') -> Dict:
    psize = np.frombuffer(fstream.read(4), np.dtype('>i4'))[0]
    flags = np.frombuffer(fstream.read(4), np.dtype('>u4'))[0]
    time, surf = np.frombuffer(fstream.read(8), np.dtype('>i4'))
    pts = np.zeros(shape=(psize, 3), dtype=np.float32) # dtype(...)
    for i in range(psize):
        pts[i, 0], pts[i, 1], pts[i, 2] = np.frombuffer(fstream.read(12), np.dtype('>f4'))
    return {
        'psize': psize, 
        'flags': flags,
        'time': time,
        'surf': surf,
        'pts': pts
    }

def parse_minx(fstream: 'BytesIO') -> Dict:
    chunk_size = np.frombuffer(fstream.read(4), np.dtype('>i4'))[0]
    oscale = np.frombuffer(fstream.read(12), np.dtype('>f4'))
    otrans = np.frombuffer(fstream.read(12), np.dtype('>f4'))
    orot = np.frombuffer(fstream.read(12), np.dtype('>f4'))
    cscale = np.frombuffer(fstream.read(12), np.dtype('>f4'))
    ctrans = np.frombuffer(fstream.read(12), np.dtype('>f4'))
    crot = np.frombuffer(fstream.read(12), np.dtype('>f4'))
    return {
        'oscale': oscale,
        'otrans': otrans,
        'orot': orot,
        'cscale': cscale,
        'ctrans': ctrans,
        'crot': crot
    }

def write_header(fstream: 'BytesIO'):
    header_data = fstream.read(236)
    return header_data

def write_objt(fstream: 'BytesIO', num_conts: int):
    name = fstream.read(64)
    extra_data = fstream.read(64)
    cont_size = np.frombuffer(fstream.read(4), np.dtype('>i4'))[0]
    new_cont_size = struct.pack('>i', num_conts)
    flags = fstream.read(4)
    axis_draw_mode = fstream.read(8)
    red_green_blue = fstream.read(12)
    pdraw_size = fstream.read(4)
    symbol_symsize_lw2D_lw3D_linesty_symflags_sympad_trans = fstream.read(8)
    mesh_size_surf_size = fstream.read(8)

    return cont_size, (name, extra_data, new_cont_size, flags, axis_draw_mode, red_green_blue, pdraw_size, symbol_symsize_lw2D_lw3D_linesty_symflags_sympad_trans, mesh_size_surf_size)

def skip_cnts(fstream: 'BytesIO'):
    psize = np.frombuffer(fstream.read(4), np.dtype('>i4'))[0]
    flags_time_surf = fstream.read(12)
    for i in range(psize):
        fstream.read(12)
    return flags_time_surf

def write_cnts(cont_num: int, contour_data: dict, flags_time_surf):
    binary_cont_data = b''
    cont_key = f'group_{cont_num}'
    psize = len(contour_data[cont_key])
    binary_cont_data = b''.join([struct.pack('>i', psize), flags_time_surf])
    for i in range(psize):
        x_bin = struct.pack('>f', contour_data[cont_key][i][0])
        y_bin = struct.pack('>f', contour_data[cont_key][i][1])
        z_bin = struct.pack('>f', contour_data[cont_key][i][2])
        binary_cont_data = b''.join([binary_cont_data, x_bin, y_bin, z_bin])
    return binary_cont_data
