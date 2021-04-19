#ifndef IMAGE_GUARD
#define IMAGE_GUARD

#include "CUDASkel2D/include/field.h"
#include "ImageWriter.hpp"
#include <set>
#include <vector>

class Image {
  private:
    /** VARIABLES **/
    float      islandThreshold;
    float             layerThreshold;
    double            *importance;
    int               numLayers;
    int               nPix; /* Short for (dimX * dimY) */
    string            compress_method;
    int*  graylevels = nullptr;
    //skel_tree_t* traceLayer(FIELD<float>* skel, FIELD<float>* dt);
    coord2D_list_t *neighbours(int x, int y, FIELD<float> *skel);

  public: 
    /** VARIABLES **/
    FIELD<float>   *im;

    /** CONSTRUCTORS **/
    Image(FIELD<float> *in);
    Image(FIELD<float> *in, float islandThresh, float importanceThresh);

    /** DESTRUCTOR **/
    ~Image();

    /** FUNCTIONS **/
    skel_tree_t* computeBinarySkel();//wang
    void encodeforcontour();//wang
    void RemoveRepeatPoints();
    void computeCUDASkeletons();
    void removeDuplicatePoints(FIELD<float> *imPrev, FIELD<float> *skP, FIELD<float> *imCur, FIELD<float> *skC);
   
};

#endif
