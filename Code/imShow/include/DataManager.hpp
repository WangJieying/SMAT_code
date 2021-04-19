/*
 * File:   DataManager.hpp
 * Author: Yuri Meiburg
 *
 * Created on May 1, 2011, 12:21 PM
 */

#ifndef DATAMANAGER_HPP
#define DATAMANAGER_HPP
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "framebuffer.hpp"
#include "texture3D.hpp"
#include "glslProgram.h"
#include "ImageDecoder.hpp"
//#include "CUDASkel2D/include/field.h"

#include <vector>
using namespace std;

typedef enum SHADER_TYPE_ {
    INVALID = -1,
    NORMAL,
    LOCAL_ALPHA,
    GLOBAL_ALPHA
} SHADER_TYPE;


typedef struct DataArr_ {
    float *data;
    int nEl;
} DataArr;

class DataManager {
private:
    int *dim;
    int clear_color;
    vector<int> gray_levels;
    Texture3D *dataLoc;
    Framebuffer *fbo;
    GLSLProgram *sh_alpha;
    SHADER_TYPE SHADER;

public:
    DataManager();
    ~DataManager();

    Texture2D* get_interp_layer(int intensity, int prev_intensity);
    void initFBO();
    void setAlphaShader(SHADER_TYPE st);
    void setData(layer_t *im) {   data = im;    }
    void setWidth (int x) { dim[0] =  x; };
    void setHeight(int y) { dim[1] =  y; };
    void set_clear_inty(int c) {clear_color = c;}
    int getWidth() {     return dim[0];    };
    int getHeight() {    return dim[1];    };
    void        set_gray_levels(vector<int>& g_l) { gray_levels = g_l;};
    vector<int>& get_gray_levels() { return gray_levels;};
    void setClearColor();
    void set_clear_color(int c) { clear_color = c;}
    FIELD<float>* get_texture_data();
    void set_fbo_data(unsigned char* texdata) {fbo->texture->setData(texdata);};
    void initCUDA();

    Texture2D* getAlphaMapOfLayer(int l);
    Texture2D* getAlphaMapOfLayer(FIELD<float>* impmap, int threshold, int superR);
private:
    //image_t *data;
    layer_t * data;
    layer_t *readLayer();
    FIELD<float>* get_alpha_to_field(int intensity);
    FIELD<float>* get_dt_of_alpha(FIELD<float>* alpha, bool background);

};

#endif  /* DATAMANAGER_HPP */

