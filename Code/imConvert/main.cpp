/* main.cpp */


#include <math.h>
#include "fileio/fileio.hpp"
#include <string>
#include <sys/resource.h>
#include "main.hpp"
#include <omp.h>
#include <vector>
#include "include/Image.hpp"
#include "include/io.hpp"
#include "include/messages.h"
#include "include/image.h"
#include "include/ImageWriter.hpp"
#include "configParser/include/Config.hpp"
#include <chrono>
#include <boost/algorithm/string.hpp>
using namespace std;


/* All properties that can be set in the config files: */
/* For a more detailed explanation, see our example.conf */
int MSG_LEVEL = MSG_NORMAL;
string filename_stdstr;
string compress_method;
/* Image space variables */
float islandThreshold = 0;
double layerThreshold = 0;
extern float hausdorff;
COLORSPACE c_space;
extern float SKELETON_DT_THRESHOLD;
extern float SKELETON_SALIENCY_THRESHOLD;
extern float IMP_THRESHOLD;
extern float SKELETON_ISLAND_THRESHOLD;
extern string OUTPUT_FILE;
extern string COMP_METHOD_NAME;
extern string ENCODING;
extern bool MSIL_SELECTION;
// Overlap pruning
extern bool OVERLAP_PRUNE;
int RemoveRepCR;
extern float OVERLAP_PRUNE_EPSILON;
// Bundling variables
extern bool BUNDLE;
extern int EPSILON;
extern float ALPHA;
extern int COMP_LEVEL;
/* Tree variables */
extern int MIN_OBJECT_SIZE;
extern int MIN_SUM_RADIUS;
extern int MIN_PATH_LENGTH;

// Set special colorspace in case of color image
// Defaults to RGB when an color image is encountered, otherwise to gray
COLORSPACE set_color_space(string c_space) {
    if (boost::iequals(c_space, "hsv")) {
        return COLORSPACE::HSV;
    } else if (boost::iequals(c_space, "ycbcr") || boost::iequals(c_space, "yuv")) {
        return COLORSPACE::YCC;
    } else if (boost::iequals(c_space, "rgb")) {
        return COLORSPACE::RGB;
    } else {
        return COLORSPACE::NONE;
    }
}
/* Set all parameters according to the configuration file parsed
 * with the Config object */
void setConfig(Config c) {
    if (!c.exists("filename")) {
        cout << "Please specify a filename" << endl;
        exit(-1);
    }
    filename_stdstr = c.pString("filename");
    if (!c.exists("outputLevel")) {
        MSG_LEVEL = MSG_NORMAL;
    } else {
        switch (c.pString("outputLevel")[0]) {
        case 'q':
            MSG_LEVEL = MSG_NEVER;
            break;
        case 'e':
            MSG_LEVEL = MSG_ERROR;
            break;
        case 'n':
            MSG_LEVEL = MSG_NORMAL;
            break;
        case 'v':
            MSG_LEVEL = MSG_VERBOSE;
            break;
        }
    }

    // MSG_LEVEL = MSG_VERBOSE;
    layerThreshold = c.exists("num_layers") ? c.pInt("num_layers") : (c.exists("lThreshold") ? c.pDouble("lThreshold") : 0);
    islandThreshold = c.exists("islandThreshold") ? c.pDouble("islandThreshold") : 0;
    COMP_METHOD_NAME = c.exists("compression_method") ? c.pString("compression_method") : "lzma";
    COMP_LEVEL = c.exists("compression_level") ? c.pInt("compression") : -1;
    ENCODING = c.exists("encoding") ? c.pString("encoding") : "standard";
    MSIL_SELECTION = c.exists("layer_selection") ? c.pString("layer_selection") != "threshold" : false;
    SKELETON_DT_THRESHOLD = c.exists("sdtThreshold") ? c.pDouble("sdtThreshold") : 5;
    SKELETON_SALIENCY_THRESHOLD = c.exists("ssThreshold") ? c.pDouble("ssThreshold") : 5;
    IMP_THRESHOLD = c.exists("Imp") ? c.pDouble("Imp") : 4;
    SKELETON_ISLAND_THRESHOLD = c.exists("SkeletonThreshold") ? c.pDouble("SkeletonThreshold") : 10;

    hausdorff = c.exists("hausdorff") ? c.pDouble("hausdorff") : 0.006;
    RemoveRepCR = c.exists("removeRepCR") ? c.pInt("removeRepCR") : 1;
    

    MIN_OBJECT_SIZE = c.exists("minObjSize") ? c.pInt("minObjSize") : 5;
    MIN_SUM_RADIUS = c.exists("minSumRadius") ? c.pInt("minSumRadius") : 15;
    MIN_PATH_LENGTH = c.exists("minPathLength") ? c.pInt("minPathLength") : 3;
 
    OVERLAP_PRUNE = c.exists("overlap_pruning") ? c.pBool("overlap_pruning") : false;
    OVERLAP_PRUNE_EPSILON = c.exists("overlap_pruning_epsilon") ? c.pDouble("overlap_pruning_epsilon") : 0.0001;
    BUNDLE = c.exists("bundle") ? c.pBool("bundle") : false;
    ALPHA = c.exists("alpha") ? min(1, max(0, c.pDouble("alpha"))) : 0;
    EPSILON = c.exists("epsilon") ? c.pInt("epsilon") : 0;
    OUTPUT_FILE = c.exists("outputFile") ? c.pString("outputFile") : "default.smat";
    c_space = c.exists("colorspace") ? set_color_space(c.pString("colorspace")) : COLORSPACE::NONE;
}


void execute_binary_pipeline(FIELD<float>* im) {
   
    float diagonal  = sqrt((float)(im->dimX() * im->dimX() + im->dimY() * im->dimY()));
    ofstream OutFile2;
    OutFile2.open("../output.txt",ios_base::app);
    OutFile2<<"importance: " << SKELETON_SALIENCY_THRESHOLD<<endl;
    OutFile2.close();
    
    int clear_color = 0;
    Image* il = new Image(im, islandThreshold, layerThreshold);
    il -> encodeforcontour();
    skel_tree_t* tree = il -> computeBinarySkel();//after this, control points have already written into the .txt file.
    
    if(RemoveRepCR) il->RemoveRepeatPoints();
    ImageWriter iw(OUTPUT_FILE.c_str());
    
    if(RemoveRepCR) iw.writeCPandID(im->dimX(), im->dimY());
    else iw.writeCP_(im->dimX(), im->dimY());
    
    delete im;
}

int main(int argc, char **argv) {
    // Increase the stacksize because default is not enough for some reason.
    const rlim_t kStackSize = 128 * 1024 * 1024;   // min stack size = 128 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0) {
        if (rl.rlim_cur < kStackSize) {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0) {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            } else {
                PRINT(MSG_NORMAL, "Set stack limit to %d\n", kStackSize);
            }
        }
    }

    /* Parse options */
    if (argc != 2) {
        fstream fhelp(". / doc / USAGE", fstream::in);
        cout << fhelp.rdbuf();
        fhelp.close();

        exit(-1);
    }
    Config config = Config(argv[1], NULL);
    setConfig(config);
    if (argc == 3) {
        OUTPUT_FILE = argv[2];
    }
    const char* filename = filename_stdstr.c_str();
    PRINT(MSG_NORMAL, "Reading input file : % s\n", filename);
    
    FIELD<float>* field = FIELD<float>::read(filename);
    if (!field) {
        PRINT(MSG_ERROR, "Failed to read file.\n");
        exit(EXIT_FAILURE);
    }
  
    execute_binary_pipeline(field);


    PRINT(MSG_NORMAL, "Done with everything.\n");
    return 0;
}




