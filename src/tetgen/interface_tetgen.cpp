#include "interface_tetgen.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include <new>

#include "../triangle/constants.h"
#include "tetgen.h"

void drop_tetgen(struct ExtTetgen *tetgen) {
    if (tetgen == NULL) {
        return;
    }
    delete tetgen;
}

struct ExtTetgen *new_tetgen(int npoint, int nfacet, int const *facet_npoint, int nregion, int nhole) {
    if (npoint < 4) {
        return NULL;
    }

    // tetgen
    struct ExtTetgen *tetgen = new (std::nothrow) ExtTetgen;
    if (tetgen == NULL) {
        return NULL;
    }
    try {
        tetgen->input.initialize();
        tetgen->output.initialize();
    } catch (...) {
        drop_tetgen(tetgen);
        return NULL;
    }

    // points
    tetgen->input.firstnumber = 0;
    tetgen->input.numberofpoints = npoint;
    tetgen->input.pointlist = new (std::nothrow) double[npoint * 3];
    if (tetgen->input.pointlist == NULL) {
        drop_tetgen(tetgen);
        return NULL;
    }

    // facets
    if (nfacet > 0) {
        tetgen->input.numberoffacets = nfacet;
        tetgen->input.facetlist = new (std::nothrow) tetgenio::facet[nfacet];
        if (tetgen->input.facetlist == NULL) {
            drop_tetgen(tetgen);
            return NULL;
        }
        const int NUM_POLY = 1;
        for (int index = 0; index < nfacet; index++) {
            // facet polygon
            tetgenio::facet *fac = &tetgen->input.facetlist[index];
            fac->polygonlist = new (std::nothrow) tetgenio::polygon[NUM_POLY];
            if (fac->polygonlist == NULL) {
                drop_tetgen(tetgen);
                return NULL;
            }
            fac->numberofpolygons = NUM_POLY;
            fac->numberofholes = 0;
            fac->holelist = NULL;
            // face polygon vertices
            size_t nvertex = facet_npoint[index];
            tetgenio::polygon *gon = &fac->polygonlist[0];
            gon->vertexlist = new (std::nothrow) int[nvertex];
            if (gon->vertexlist == NULL) {
                drop_tetgen(tetgen);
                return NULL;
            }
            gon->numberofvertices = nvertex;
        }
    }

    // regions
    if (nregion > 0) {
        tetgen->input.numberofregions = nregion;
        tetgen->input.regionlist = new (std::nothrow) double[nregion * 5];
        if (tetgen->input.regionlist == NULL) {
            drop_tetgen(tetgen);
            return NULL;
        }
    }

    // holes
    if (nhole > 0) {
        tetgen->input.numberofholes = nhole;
        tetgen->input.holelist = new (std::nothrow) double[nhole * 3];
        if (tetgen->input.holelist == NULL) {
            drop_tetgen(tetgen);
            return NULL;
        }
    }

    return tetgen;
}

int tet_set_point(struct ExtTetgen *tetgen, int index, double x, double y, double z) {
    if (tetgen == NULL) {
        return TRITET_ERROR_NULL_DATA;
    }
    if (tetgen->input.pointlist == NULL) {
        return TRITET_ERROR_NULL_POINT_LIST;
    }
    if (index >= tetgen->input.numberofpoints) {
        return TRITET_ERROR_INVALID_POINT_INDEX;
    }
    tetgen->input.pointlist[index * 3] = x;
    tetgen->input.pointlist[index * 3 + 1] = y;
    tetgen->input.pointlist[index * 3 + 2] = z;

    return TRITET_SUCCESS;
}

int tet_set_facet_point(struct ExtTetgen *tetgen, int index, int m, int p) {
    if (tetgen == NULL) {
        return TRITET_ERROR_NULL_DATA;
    }
    if (tetgen->input.facetlist == NULL) {
        return TRITET_ERROR_NULL_FACET_LIST;
    }
    if (index >= tetgen->input.numberoffacets) {
        return TRITET_ERROR_INVALID_FACET_INDEX;
    }

    tetgenio::facet *fac = &tetgen->input.facetlist[index];
    if (fac->polygonlist == NULL) {
        return TRITET_ERROR_NULL_FACET_POLYGON_LIST;
    }
    if (fac->numberofpolygons != 1) {
        return TRITET_ERROR_INVALID_FACET_NUM_POLYGON;
    }

    tetgenio::polygon *gon = &fac->polygonlist[0];
    if (m >= gon->numberofvertices) {
        return TRITET_ERROR_INVALID_FACET_POINT_INDEX;
    }
    if (p >= tetgen->input.numberofpoints) {
        return TRITET_ERROR_INVALID_FACET_POINT_ID;
    }
    gon->vertexlist[m] = p;

    return TRITET_SUCCESS;
}

int tet_set_region(struct ExtTetgen *tetgen, int index, double x, double y, double z, int attribute, double max_volume) {
    if (tetgen == NULL) {
        return TRITET_ERROR_NULL_DATA;
    }
    if (tetgen->input.regionlist == NULL) {
        return TRITET_ERROR_NULL_REGION_LIST;
    }
    if (index >= tetgen->input.numberofregions) {
        return TRITET_ERROR_INVALID_REGION_INDEX;
    }
    tetgen->input.regionlist[index * 5] = x;
    tetgen->input.regionlist[index * 5 + 1] = y;
    tetgen->input.regionlist[index * 5 + 2] = z;
    tetgen->input.regionlist[index * 5 + 3] = attribute;
    tetgen->input.regionlist[index * 5 + 4] = max_volume;

    return TRITET_SUCCESS;
}

int tet_set_hole(struct ExtTetgen *tetgen, int index, double x, double y, double z) {
    if (tetgen == NULL) {
        return TRITET_ERROR_NULL_DATA;
    }
    if (tetgen->input.holelist == NULL) {
        return TRITET_ERROR_NULL_HOLE_LIST;
    }
    if (index >= tetgen->input.numberofholes) {
        return TRITET_ERROR_INVALID_HOLE_INDEX;
    }
    tetgen->input.holelist[index * 3] = x;
    tetgen->input.holelist[index * 3 + 1] = y;
    tetgen->input.holelist[index * 3 + 2] = z;

    return TRITET_SUCCESS;
}

int tet_run_delaunay(struct ExtTetgen *tetgen, int verbose) {
    if (tetgen == NULL) {
        return TRITET_ERROR_NULL_DATA;
    }
    if (tetgen->input.pointlist == NULL) {
        return TRITET_ERROR_NULL_POINT_LIST;
    }

    // Tetrahedralize the points
    // Switches:
    // * `z` -- number everything from zero (z)
    char command[10];
    strcpy(command, "z");
    if (verbose == TRITET_FALSE) {
        strcat(command, "Q");
    }
    try {
        tetrahedralize(command, &tetgen->input, &tetgen->output, NULL, NULL);
    } catch (int status) {
        printf("status = %d\n", status);  // TODO
    } catch (...) {
        return 1;  // TODO
    }

    return TRITET_SUCCESS;
}

int tet_run_tetrahedralize(struct ExtTetgen *tetgen, int verbose, int o2, double global_max_volume, double global_min_angle) {
    if (tetgen == NULL) {
        return TRITET_ERROR_NULL_DATA;
    }
    if (tetgen->input.pointlist == NULL) {
        return TRITET_ERROR_NULL_POINT_LIST;
    }
    if (tetgen->input.facetlist == NULL) {
        return TRITET_ERROR_NULL_FACET_LIST;
    }

    // Generate mesh
    // Switches:
    // * `p` -- tetrahedralize a piecewise linear complex (PLC)
    // * `z` -- number everything from zero (z)
    // * `A` -- assign a regional attribute to each element (A)
    char command[128];
    strcpy(command, "pzA");
    if (verbose == TRITET_FALSE) {
        strcat(command, "Q");
    }
    if (o2 == TRITET_TRUE) {
        strcat(command, "o2");
    }
    if (global_max_volume > 0.0) {
        char buf[32];
        int n = snprintf(buf, 32, "a%.15f", global_max_volume);
        if (n >= 32) {
            return TRITET_ERROR_STRING_CONCAT;
        }
        strcat(command, buf);
    }
    if (global_min_angle > 0.0) {
        char buf[32];
        int n = snprintf(buf, 32, "q%.15f", global_min_angle);
        if (n >= 32) {
            return TRITET_ERROR_STRING_CONCAT;
        }
        strcat(command, buf);
    } else {
        strcat(command, "q");
    }
    try {
        tetrahedralize(command, &tetgen->input, &tetgen->output, NULL, NULL);
    } catch (int status) {
        printf("status = %d\n", status);  // TODO
    } catch (...) {
        return 1;  // TODO
    }

    return TRITET_SUCCESS;
}

int tet_get_npoint(struct ExtTetgen *tetgen) {
    if (tetgen == NULL) {
        return 0;
    }
    return tetgen->output.numberofpoints;

    return 0;
}

int tet_get_ntetrahedron(struct ExtTetgen *tetgen) {
    if (tetgen == NULL) {
        return 0;
    }
    return tetgen->output.numberoftetrahedra;

    return 0;
}

int tet_get_ncorner(struct ExtTetgen *tetgen) {
    if (tetgen == NULL) {
        return 0;
    }
    return tetgen->output.numberofcorners;

    return 0;
}

double tet_get_point(struct ExtTetgen *tetgen, int index, int dim) {
    if (tetgen == NULL) {
        return 0.0;
    }
    if (index < tetgen->output.numberofpoints && (dim == 0 || dim == 1 || dim == 2)) {
        return tetgen->output.pointlist[index * 3 + dim];
    } else {
        return 0.0;
    }

    return 0.0;
}

int tet_get_tetrahedron_corner(struct ExtTetgen *tetgen, int index, int corner) {
    if (tetgen == NULL) {
        return 0;
    }
    if (index < tetgen->output.numberoftetrahedra && corner < tetgen->output.numberofcorners) {
        return tetgen->output.tetrahedronlist[index * tetgen->output.numberofcorners + corner];
    } else {
        return 0;
    }

    return 0;
}

int tet_get_tetrahedron_attribute(struct ExtTetgen *tetgen, int index) {
    if (tetgen == NULL) {
        return 0;
    }
    if (index < tetgen->output.numberoftetrahedra && tetgen->output.numberoftetrahedronattributes > 0) {
        return tetgen->output.tetrahedronattributelist[index * tetgen->output.numberoftetrahedronattributes];
    } else {
        return 0;
    }

    return 0;
}
