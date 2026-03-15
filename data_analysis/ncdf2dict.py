"""
ncdf2dict.py
-------------
Robust NetCDF to Python dictionary reader, with support for GS2-style complex variables.
"""

import numpy as np

try:
    from netCDF4 import Dataset
except ImportError:
    raise ImportError("ncdf2dict.py requires netCDF4 module. Please install it via pip.")

def ncdf2dict(filename, comDim='ri', only=None, include_global_attrs=True):
    """
    Read a NetCDF file and convert it to a Python dictionary.

    Parameters
    ----------
    filename : str
        Path to the NetCDF file.
    comDim : str, optional
        Name of the dimension representing the real/imag parts of complex numbers (default 'ri').
    only : list of str, optional
        If provided, only include these variable names.
    include_global_attrs : bool, optional
        If True, include global attributes in the dictionary under 'global_attributes'.

    Returns
    -------
    dict
        Dictionary mapping variable names to NumPy arrays, with complex numbers reconstructed if comDim is present.
    """
    if filename is None:
        raise ValueError("Must provide a NetCDF filename")

    data_dict = {}
    ds = None

    try:
        ds = Dataset(filename, 'r')
        for varname in ds.variables:
            if only is not None and varname not in only:
                continue
            var = ds.variables[varname]
            dims = var.dimensions
            # Check if this variable has a complex dimension
            if comDim in dims:
                dim_index = dims.index(comDim)
                if var.shape[dim_index] != 2:
                    raise ValueError(f"Expected complex dimension '{comDim}' to have size 2 in variable '{varname}'")
                # Move complex dimension to last
                order = list(range(len(dims)))
                order.append(order.pop(dim_index))
                data = var[:].transpose(order)
                shape_new = list(data.shape)
                shape_new.pop(-1)
                # Combine real+imag
                data_flat = data.reshape(-1, 2)
                data_complex = data_flat[:,0] + 1j * data_flat[:,1]
                data_dict[varname] = data_complex.reshape(shape_new)
            else:
                data_dict[varname] = var[:]

        if include_global_attrs:
            attrs = {name: getattr(ds, name) for name in ds.ncattrs()}
            data_dict['global_attributes'] = attrs

    except Exception as e:
        print(f"ERROR reading {filename}: {e}")
        return None

    finally:
        if ds is not None:
            ds.close()

    return data_dict

