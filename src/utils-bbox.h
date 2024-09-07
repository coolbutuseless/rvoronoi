

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Bounding box type
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
typedef struct {
  double xmin;
  double ymin;
  double xmax;
  double ymax;
} bbox_t;

bbox_t bbox_new(void);
void bbox_add(bbox_t *bbox, int n, double *x, double *y);
void bbox_expand(bbox_t *bbox, double factor);
