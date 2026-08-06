"""Small grid utilities for the dye-tracer figure scripts.

`create_coordinate_edges` reproduces the identically-named helper from the
`mw_protocol` toolbox (Rome et al., 2022), used here to build cell-edge
coordinates for `pcolormesh` plotting. It is a generic, two-line grid utility;
it is reproduced locally (rather than importing the full external package) so
these figure scripts run self-contained. The scientific meltwater-routing
methods of `mw_protocol` are NOT reproduced here — see the README for the
external dependency and its citation:

    Olnavy (2022). Olnavy/ROME2022_paleoceanography_oscillations:
    Reviews round 1 - v1.4 [Software]. Zenodo.
    https://doi.org/10.5281/zenodo.6788389
"""


def create_coordinate_edges(coordinates):
    """Return cell-edge coordinates for a 1-D, regularly spaced axis.

    Given N regularly spaced cell-centre coordinates, return the N+1 cell
    edges (each edge half a step outside the outermost centres). Suitable for
    passing as the x/y edge arrays to ``matplotlib.pyplot.pcolormesh``.

    Parameters
    ----------
    coordinates : sequence of float
        1-D, regularly spaced coordinate centres (e.g. longitudes or latitudes).

    Returns
    -------
    list of float
        N+1 edge coordinates.
    """
    step = coordinates[1] - coordinates[0]
    return [coordinates[0] - step / 2 + i * step for i in range(len(coordinates) + 1)]
