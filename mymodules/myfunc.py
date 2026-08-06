# ------------------------------------- #
# Laura adds here more useful functions and classes


# ------------------------------------- #
# ---------- GENERAL METHODS ---------- #
# ------------------------------------- #


def re_round(li, prec=5):
    """
    For a list of touples, round each element
    :param list
    :param prec: gives the precision
    :return:

    """
    try:
         return round(li, prec)
    except TypeError:
         return type(li)(re_round(x, prec) for x in li)

import time
class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print('[%s]' % self.name,)
        print('Elapsed: %s' % (time.time() - self.tstart))


### Calculate Surface area for grids
import numpy as np
def cell_area(n_lon, lat1, lat2):
    """
    Area of a cell on a regular lon-lat grid.
    :param n_lon: number of longitude divisions
    :param lat1: bottom of the cell
    :param lat2: top of the cell
    :return:
    """
    r = 6371000
    lat1_rad, lat2_rad = 2 * np.pi * lat1 / 360, 2 * np.pi * lat2 / 360
    return 2 * np.pi * r ** 2 * np.abs(np.sin(lat1_rad) - np.sin(lat2_rad)) / n_lon


def cell_volume(lon1,lon2, lat1, lat2,d1,d2):
    """
    Volume of a cell on a regular lon-lat grid (Formula Kelly & Savric, 2020).
    Result is given in m^3 (Convert to km^3 by dividing with 10^9).
    :param lon 1: lower edge of lon segment
    :param lon 2: higher edge of lon segment
    :param lat1: bottom of the cell
    :param lat2: top of the cell
    :param d1: upper depth level (in m)
    :param d2: lower depth level (in m)
    :return:
    """
    r1 = 6371000 - d1
    r2 = 6371000 - d2
    lat1_rad, lat2_rad = 2 * np.pi * lat1 / 360, 2 * np.pi * lat2 / 360
    lon1_rad, lon2_rad = 2 * np.pi * lon1 / 360, 2 * np.pi * lon2 / 360

    return 1/3*(np.abs((r2**3-r1**3))*np.abs(lon2_rad - lon1_rad)*np.abs(np.sin(lat2_rad)-np.sin(lat1_rad))) 



def guess_bounds(coordinate):
    if coordinate is not None:
        if len(coordinate) <= 1:
            coordinateb = coordinate
        else:
            coordinateb = [(coordinate[i] + coordinate[i + 1]) / 2 for i in range(len(coordinate) - 1)]
            coordinateb = np.append((3 * coordinate[0] - coordinate[1]) / 2, coordinateb)
            coordinateb = np.append(coordinateb, (3 * coordinate[-1] - coordinate[-2]) / 2)
        return np.array(coordinateb)
    else:
        raise ValueError("Empty coordinate.")

def surface_matrix(lon, lat):
    """
    Compute a matrix with all the surfaces values.
    :param lon:
    :param lat:
    :return:
    """
    n_j, n_i = len(lat), len(lon)
    lat_b = guess_bounds(lat)
    surface = np.zeros((n_j, n_i))
    for i in range(n_i):
        for j in range(n_j):
            surface[j, i] = cell_area(n_i, lat_b[j], lat_b[j + 1])
    return surface


### -----------------
### Related to d18O
### -----------------



def build_d18O_lookup(df, column, dye_index=None):
    """
    Build dye → d18O lookup dictionary from a pandas DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Your gdg dataframe

    column : str
        Column containing d18O values (or region/scenario column)

    dye_index : list or None
        Optional list of dye names (dye00...dye08).
        If None, assumes 9 dyes in order.
    """

    import numpy as np
    import xarray as xr

    if dye_index is None:
        dye_index = [f"dye{i:02d}" for i in range(9)]

    values = df[column].values

    if len(values) != len(dye_index):
        raise ValueError(
            f"Mismatch: {len(values)} values vs {len(dye_index)} dyes"
        )

    return dict(zip(dye_index, values))

import numpy as np
import xarray as xr


def compute_d18O(ds, lookup, dyes=None):
    """
    Compute total d18O field from dye tracers.

    Parameters
    ----------
    ds : xarray.Dataset
        Must contain dye variables: dye00 ... dye08 (or subset)

    lookup : dict
        Mapping dye -> d18O value

    dyes : list (optional)
        Explicit dye order. If None, inferred from lookup keys.

    Returns
    -------
    xarray.DataArray
        Total d18O field (same spatial/time dimensions as dyes)
    """

    if dyes is None:
        dyes = list(lookup.keys())

    # --- stack dyes into a single dimension
    dye_stack = xr.concat(
        [ds[d] for d in dyes],
        dim=xr.IndexVariable("dye", dyes)
    )

    # --- weights (d18O per dye)
    weights = xr.DataArray(
        np.array([lookup[d] for d in dyes]),
        dims="dye",
        coords={"dye": dyes},
    )

    # --- weighted sum
    return (dye_stack * weights).sum("dye")


