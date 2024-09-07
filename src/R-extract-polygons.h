SEXP extract_polygons_internal(
    int vert_n, double *vert_x, double *vert_y, // voronoi vertices
    int seg_n, int *seg_v1, int *seg_v2,        // voronoi edges
    int seed_n, double *seed_x, double *seed_y  // seed points (if sorting is required)
);
