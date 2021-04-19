#include "include/DataManager.hpp"
#include "include/texture3D.hpp"
#include "include/framebuffer.hpp"
#include "include/io.hpp"
#include "include/messages.h"
#include "include/skeleton_cuda.hpp"
#include <math.h>
#define ARRAY_SIZE(arr) sizeof(arr)/sizeof(arr[0])
#define SET_TEXTURE(arr, i, val) do { (arr)[(4 * (i))] = (arr)[(4 * (i) + 1)] = (arr)[(4 * (i) + 2)] = (val); (arr)[(4 * (i) + 3)] = 255.0; } while(0);
DataManager::DataManager() {
    dim = new int[2];
    fbo = (Framebuffer *) 0;
}

DataManager::~DataManager() {
    delete [] dim;
    delete dataLoc;
    delete fbo;
    delete sh_alpha;
    deallocateCudaMem();
}

void DataManager::initCUDA() {
    initialize_skeletonization(dim[0], dim[1]);

}

// Blit the texture to a field so we can work with it.
FIELD<float>* DataManager::get_alpha_to_field(int intensity) {
    if (intensity == 0)
        return nullptr;
    float *data = (float *) malloc(dim[0] * dim[1] * sizeof (float));
    Texture2D *tex = getAlphaMapOfLayer(intensity);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, tex->tex);

    /* Altering range [0..1] -> [0 .. 255] */
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, data);
 
    FIELD<float> *im = new FIELD<float>(0, 0);
    im->setAll(dim[0], dim[1], data);
   
    return im;
}

FIELD<float>* DataManager::get_texture_data() {
    unsigned char* data = (unsigned char*)(calloc(dim[0] * dim[1], sizeof(unsigned char)));
    glReadPixels(0, 0, dim[0], dim[1], GL_RED, GL_UNSIGNED_BYTE, data);
    FIELD<float>* f = new FIELD<float>(dim[0], dim[1]);
    for (int i = 0; i < dim[0] * dim[1]; ++i) {
        int y = i / dim[0];
        int x = i % dim[0];
        f->set(x, y, data[i]);
    }
    return f;
}

// Computes the DT of the alpha map. It's settable if the DT needs to be computed
// of the foreground (i.e. pixels which are 0) or the background (i.e. pixels that are 1)
FIELD<float>* DataManager::get_dt_of_alpha(FIELD<float>* alpha, bool foreground) {
    auto alpha_dupe = alpha->dupe();
    auto data = alpha_dupe->data();
    // Widen the image range from [0..1] to [0..255]

    for (int i = 0; i < dim[0] * dim[1]; i++) {
        data[i] *= 255;
    }
    
    auto ret = computeCUDADT(alpha_dupe, foreground);
   
    return ret;
}

Texture2D* DataManager::get_interp_layer(int intensity, int prev_intensity) {
    
    unsigned char *texdata = (unsigned char *) calloc(dim[0] * dim[1] * 4, sizeof (char));
    bool compute_dt_of_foreground = false;

    // Compute all data
    FIELD<float>* curr_alpha = get_alpha_to_field(intensity);
    
    FIELD<float>* curr_dt = get_dt_of_alpha(curr_alpha, !compute_dt_of_foreground);

    curr_alpha->writePGM(string("temp/" + to_string(intensity) + "a.pgm").c_str());
    curr_dt->writePGM(string("temp/" + to_string(intensity) + "b.pgm").c_str());

    FIELD<float>* prev_alpha = nullptr;
    FIELD<float>* prev_dt = nullptr;
    float* prev_alpha_data = nullptr;
    float* prev_dt_data = nullptr;
    float* curr_alpha_data = curr_alpha->data();
    float* curr_dt_data = curr_dt->data();

    // If there is a previous layer we might be able to interpolate
    if (prev_intensity != 0) {
        // Alpha layer of previous threshold and DT of its foreground
        prev_alpha = get_alpha_to_field(prev_intensity);
        prev_dt = get_dt_of_alpha(prev_alpha, compute_dt_of_foreground);

        prev_alpha->writePGM(string("temp/" + to_string(intensity) + "c.pgm").c_str());
        prev_dt->writePGM(string("temp/" + to_string(intensity) + "d.pgm").c_str());
        prev_alpha_data = prev_alpha->data();
        prev_dt_data = prev_dt->data();

        for (int i = 0; i < dim[0] * dim[1]; ++i) {
            // If the current foreground is set we set it to that
            if (curr_alpha_data[i]) {
                SET_TEXTURE(texdata, i, intensity);
            } else {
                // If there are pixels active between boundaries we smoothly interpolate them
                if (prev_alpha_data[i]) {
                    float prev_dt_val = prev_dt_data[i];
                    float curr_dt_val = curr_dt_data[i];

                    float interp_color = 0.5 * (min(1, prev_dt_val / curr_dt_val) * prev_intensity + max(1 - curr_dt_val / prev_dt_val, 0) * intensity);
                    SET_TEXTURE(texdata, i, interp_color);
                }
            }
            // Otherwise we keep the previous set value
        }
    } else {
        // Set everything to current layer
        for (int i = 0; i < dim[0] * dim[1]; ++i) {
            SET_TEXTURE(texdata, i, intensity);
        }
    }
    delete curr_alpha;
    delete prev_alpha;
    delete curr_dt;
    delete prev_dt;
    fbo->texture->setData(texdata);
    free(texdata);
    return fbo->texture;
}
/* This function draws all disks. The alpha will be 1 at the location that should be drawn, otherwise 0. */
Texture2D* DataManager::getAlphaMapOfLayer(FIELD<float>* impmap,int threshold, int superR) {
    
    //impmap->writePGM("impppp.pgm");
    FIELD<float>* skel = new FIELD<float>(dim[0], dim[1]);
    for (unsigned int x = 0; x < dim[0]; ++x) 
        for (unsigned int y = 0; y < dim[1]; ++y)
            skel->set(x, y, 0);
    

    layer_t *da = readLayer();
    int x, y, radius;

    if (fbo == 0) {
        PRINT(MSG_ERROR, "Framebuffer object was not initialized. Creating new FBO (FIX THIS!)\n");
        initFBO();
    }

    /* Draw on a texture. */
    fbo->bind();
    if (SHADER == GLOBAL_ALPHA) {
        glClearColor(1.0, 1.0, 1.0, 1.0);
    }

    glEnable(GL_DEPTH_TEST);
    glClear(GL_DEPTH_BUFFER_BIT);
    glClear(GL_COLOR_BUFFER_BIT);
    sh_alpha->bind();

    glBegin(GL_QUADS);
    for (unsigned int j = 0; j < da->size(); ++j) {
        x = (*da)[j].first;
        if( dim[0]>1920 || dim[1]>1080) y = dim[1] - (*da)[j].second;
        else y = (*da)[j].second;
        radius = (*da)[j].third;

        //if (impmap->value(x,y)<threshold) continue;//4D
        skel->set(x, y, 1);
         // Texture coordinates are in range [-1, 1]
        // Texcoord selects the texture coordinate
        // vertex call sets it to that location
     
    //x = 100; y=100;  radius = 2;
        glTexCoord2f(-1.0, -1.0);
        glVertex2f(x - radius, y - radius);
        glTexCoord2f(-1.0, 1.0);
        glVertex2f(x - radius, y + radius+1);
        glTexCoord2f(1.0, 1.0);
        glVertex2f(x + radius+1, y + radius+1);
        glTexCoord2f(1.0, -1.0);
        glVertex2f(x + radius+1, y - radius);
        
       /*
        glTexCoord2f(-1.0, -1.0);
        glVertex2f(x - radius, y - radius);
        glTexCoord2f(-1.0, 1.0);
        glVertex2f(x - radius, y + radius);
        glTexCoord2f(1.0, 1.0);
        glVertex2f(x + radius, y + radius);
        glTexCoord2f(1.0, -1.0);
        glVertex2f(x + radius, y - radius);
        */
    }
    glEnd();


    //skel->set(0, 0, 0);//I don't know why there is a skeleton point in the left top(x=0,y=0) corner, but just delete it.   
    skel->writePGM("skel.pgm");
    sh_alpha->unbind();
    fbo->unbind();

    glClearColor(1.0, 1.0, 1.0, 1.0);
    /////////sub-pixel
    if (superR>1 || dim[0]>1920 || dim[1]>1080)
    {
        float *data = (float *) malloc(dim[0] * dim[1] * sizeof (float));
        //Texture2D *tex = getAlphaMapOfLayer(intensity);

        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, fbo->texture->tex);

        // Altering range [0..1] -> [0 .. 255] 
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, data);
    
        FIELD<float> *im = new FIELD<float>(0, 0);
        im->setAll(dim[0], dim[1], data);
        im->writePGM("buffer.pgm");
    }
    return da->size() != 0 ? fbo->texture : NULL;
}
/* This function draws all disks. The alpha will be 1 at the location that should be drawn, otherwise 0. */
Texture2D* DataManager::getAlphaMapOfLayer(int l) {
    layer_t *da = readLayer();
    int x, y, radius;

    if (fbo == 0) {
        PRINT(MSG_ERROR, "Framebuffer object was not initialized. Creating new FBO (FIX THIS!)\n");
        initFBO();
    }

    /* Draw on a texture. */
    fbo->bind();
    if (SHADER == GLOBAL_ALPHA) {
        glClearColor(1.0, 1.0, 1.0, 1.0);
    }

    glEnable(GL_DEPTH_TEST);
    glClear(GL_DEPTH_BUFFER_BIT);
    glClear(GL_COLOR_BUFFER_BIT);
    sh_alpha->bind();
    glBegin(GL_QUADS);
    for (unsigned int j = 0; j < da->size(); ++j) {
        x = (*da)[j].first;
        y = (*da)[j].second;
        radius = (*da)[j].third;

        // Texture coordinates are in range [-1, 1]
        // Texcoord selects the texture coordinate
        // vertex call sets it to that location
        glTexCoord2f(-1.0, -1.0);
        glVertex2f(x - radius, y - radius);
        glTexCoord2f(-1.0, 1.0);
        glVertex2f(x - radius, y + radius);
        glTexCoord2f(1.0, 1.0);
        glVertex2f(x + radius, y + radius);
        glTexCoord2f(1.0, -1.0);
        glVertex2f(x + radius, y - radius);
    }
    glEnd();

    sh_alpha->unbind();
    fbo->unbind();

    glClearColor(0.0, 0.0, 0.0, 1.0);
    return da->size() != 0 ? fbo->texture : NULL;
}

/* Return all disks for a certain intensity */
layer_t *DataManager::readLayer() {
    return data;
}

/* When the lowest intensity is >0, we should set the clear color to i-1, where
* i is the first intensity with disks.
*/
void DataManager::setClearColor() {
   
    glClearColor(clear_color / 255.0, clear_color / 255.0, clear_color / 255.0, 0);
   
    PRINT(MSG_VERBOSE, "Clear color set to: %f, %f, %f (layer %d)\n", clear_color / 255.0, clear_color / 255.0, clear_color / 255.0, clear_color);
}

void DataManager::initFBO() {
    //fbo = new Framebuffer(1920, 1080, GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, true);
    fbo = new Framebuffer(dim[0], dim[1], GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, true);
}

void DataManager::setAlphaShader(SHADER_TYPE st) {
    SHADER = st;

    if (st == LOCAL_ALPHA) {
        sh_alpha = new GLSLProgram("glsl/alpha.vert", "glsl/alpha.frag");
        sh_alpha->compileAttachLink();
        PRINT(MSG_VERBOSE, "Succesfully compiled shader 1 (local alpha computation)\n");
    } else if (st == GLOBAL_ALPHA) {
        sh_alpha = new GLSLProgram("glsl/nointerpolationinv.vert", "glsl/nointerpolationinv.frag");
        sh_alpha->compileAttachLink();
        PRINT(MSG_VERBOSE, "Succesfully compiled shader 1 (global alpha computation)\n");
    } else if (st == NORMAL) {
        sh_alpha = new GLSLProgram("glsl/nointerpolation.vert", "glsl/nointerpolation.frag");
        sh_alpha->compileAttachLink();
        PRINT(MSG_VERBOSE, "Succesfully compiled shader 1 (No interpolation)\n");
    } else {
        PRINT(MSG_ERROR, "Invalid shader. Did not compile.");
        exit(1);
    }
}