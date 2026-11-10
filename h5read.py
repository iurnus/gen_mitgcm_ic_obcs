"""Minimal HDF5 (NetCDF-4) reader via ctypes on the system libhdf5.
Used because no Python NetCDF/HDF5 package is installed. If you have xarray/netCDF4,
the reading part of make_ic_obc.py can be swapped for them."""
import ctypes as C, ctypes.util, numpy as np

def _load():
    for name in ['libhdf5_serial.so.103', 'libhdf5.so.103', ctypes.util.find_library('hdf5_serial'),
                 ctypes.util.find_library('hdf5')]:
        if name:
            try: return C.CDLL(name)
            except OSError: pass
    raise OSError('libhdf5 not found')
lib = _load()
lib.H5open()
hid = C.c_int64
for f in ['H5Fopen', 'H5Dopen2', 'H5Dget_space', 'H5Dget_type', 'H5Aopen', 'H5Aget_type',
          'H5Aget_space', 'H5Tget_native_type']:
    getattr(lib, f).restype = hid
lib.H5Fopen.argtypes = [C.c_char_p, C.c_uint, hid]
lib.H5Dopen2.argtypes = [hid, C.c_char_p, hid]
lib.H5Dget_space.argtypes = [hid]; lib.H5Dget_type.argtypes = [hid]
lib.H5Sget_simple_extent_ndims.argtypes = [hid]
lib.H5Sget_simple_extent_dims.argtypes = [hid, C.c_void_p, C.c_void_p]
lib.H5Tget_native_type.argtypes = [hid, C.c_int]
lib.H5Tget_class.argtypes = [hid]; lib.H5Tget_size.argtypes = [hid]; lib.H5Tget_size.restype = C.c_size_t
lib.H5Tget_sign.argtypes = [hid]; lib.H5Tis_variable_str.argtypes = [hid]
lib.H5Dread.argtypes = [hid, hid, hid, hid, hid, C.c_void_p]
lib.H5Aopen.argtypes = [hid, C.c_char_p, hid]
lib.H5Aget_type.argtypes = [hid]; lib.H5Aget_space.argtypes = [hid]
lib.H5Aread.argtypes = [hid, hid, C.c_void_p]
lib.H5Aexists.argtypes = [hid, C.c_char_p]
lib.H5Lexists.argtypes = [hid, C.c_char_p, hid]
for f in ['H5Dclose', 'H5Sclose', 'H5Tclose', 'H5Aclose', 'H5Fclose']:
    getattr(lib, f).argtypes = [hid]

def _dtype(t):
    cls = lib.H5Tget_class(t); sz = lib.H5Tget_size(t)
    if cls == 1: return np.dtype(f'f{sz}')
    if cls == 0: return np.dtype(('i' if lib.H5Tget_sign(t) == 1 else 'u') + str(sz))
    raise TypeError(f'unsupported class {cls}')

def _shape(s):
    n = lib.H5Sget_simple_extent_ndims(s)
    d = (C.c_uint64 * max(n, 1))(); lib.H5Sget_simple_extent_dims(s, d, None)
    return tuple(d[i] for i in range(n))

class H5File:
    def __init__(self, path):
        self.path = path
        self.f = lib.H5Fopen(path.encode(), 0, 0)
        if self.f < 0: raise IOError(path)

    def attr(self, obj, name):
        if lib.H5Aexists(obj, name.encode()) <= 0: return None
        a = lib.H5Aopen(obj, name.encode(), 0); t = lib.H5Aget_type(a); s = lib.H5Aget_space(a)
        try:
            if lib.H5Tget_class(t) == 3:                       # string
                if lib.H5Tis_variable_str(t) > 0:
                    buf = (C.c_char_p * 1)(); lib.H5Aread(a, t, buf)
                    return buf[0].decode(errors='replace') if buf[0] else ''
                n = lib.H5Tget_size(t); buf = C.create_string_buffer(n)
                lib.H5Aread(a, t, buf); return buf.raw.split(b'\0')[0].decode(errors='replace')
            nt = lib.H5Tget_native_type(t, 0)
            try:
                arr = np.empty(_shape(s) or (1,), _dtype(nt))
                lib.H5Aread(a, nt, arr.ctypes.data)
            finally:
                lib.H5Tclose(nt)
            return arr.ravel()[0] if arr.size == 1 else arr
        except TypeError:
            return None
        finally:
            lib.H5Tclose(t); lib.H5Sclose(s); lib.H5Aclose(a)

    def read(self, name, decode=True, dtype='f8'):
        """Return (array, attrs). Applies _FillValue -> NaN, scale_factor, add_offset."""
        d = lib.H5Dopen2(self.f, name.encode(), 0)
        if d < 0: raise KeyError(f'{name} not in {self.path}')
        t = lib.H5Dget_type(d); nt = lib.H5Tget_native_type(t, 0); dt = _dtype(nt)
        s = lib.H5Dget_space(d); arr = np.empty(_shape(s), dt)
        ok = lib.H5Dread(d, nt, 0, 0, 0, arr.ctypes.data)
        at = {k: self.attr(d, k) for k in ['_FillValue', 'missing_value', 'scale_factor',
                                            'add_offset', 'units']}
        lib.H5Tclose(nt); lib.H5Tclose(t); lib.H5Sclose(s); lib.H5Dclose(d)
        if ok < 0: raise IOError('read failed ' + name)
        if decode:
            raw = arr; arr = arr.astype(dtype)
            for k in ('_FillValue', 'missing_value'):
                if at[k] is not None and not isinstance(at[k], str):
                    arr[raw == at[k]] = np.nan
            if at['scale_factor'] is not None: arr *= np.dtype(dtype).type(at['scale_factor'])
            if at['add_offset'] is not None: arr += np.dtype(dtype).type(at['add_offset'])
        return arr, at

    def has(self, name):
        return lib.H5Lexists(self.f, name.encode(), 0) > 0

    def close(self):
        lib.H5Fclose(self.f)
