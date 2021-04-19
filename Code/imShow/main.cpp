#include <fstream>
#include <cmath>
#include <sys/stat.h>

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <unistd.h>
#include <tuple>
#include "include/main.hpp"
#include "include/io.hpp"
#include "include/messages.h"
#include "include/glslProgram.h"
#include "include/framebuffer.hpp"
#include "include/DataManager.hpp"
#include "include/ImageDecoder.hpp"
#include "lodepng/lodepng.h"
#include "include/image.h"
#include "fileio/fileio.hpp"
//wang.

#include "SplineGenerate/BSplineCurveFitter/BSplineCurveFitterWindow3.h"
#include <sstream>
#include <string.h>
using namespace std;

#define MIN(a,b) (a) < (b) ? (a):(b)
#define MAX(a,b) (a) > (b) ? (a):(b)
#define ARRAY_SIZE(arr) sizeof(arr)/sizeof(arr[0])


SHADER_TYPE SHADER = NORMAL;

char * i; bool SaveMore = 0; FIELD<float>* impmap; int SuperResolution = 1;//wang
int MSG_LEVEL = MSG_VERBOSE;
int WWIDTH = 0, WHEIGHT = 0;
int clear_color;
unsigned char* new_pixel_data = 0;
DataManager *dm, *dm_g, *dm_b;
bool is_color_image = false;
COLORSPACE c_space = COLORSPACE::NONE;
GLSLProgram *sh_render;
bool interpolate = true;
int display1;
using namespace std;

char *outFile = 0; /* Press 's' to make a screenshot to location specified by outFile. */

void saveOutput() {
    unsigned char *sdata = (unsigned char *) malloc(WWIDTH * WHEIGHT * 4);
    glReadPixels(0, 0, WWIDTH, WHEIGHT, GL_RGB, GL_UNSIGNED_BYTE, sdata);
    
    Texture2D *t = new Texture2D(WWIDTH, WHEIGHT, GL_RGB, GL_RGB, GL_UNSIGNED_BYTE);
    t->setData(sdata);
    t->saveAsPNG(outFile);
    free(sdata);
    delete t;
}

void initDataManager(const char *file) {
    string str;
    ifstream ifs("config.txt");
    ifs >> str;
    ifs >> str;
    int removeRepCP = (int)atof(str.c_str());

    ImageDecoder id;
    if(removeRepCP) id.loadSMAT_R(file);
    else id.loadSMAT_(file);
      
    BSplineCurveFitterWindow3 spline;
    spline.SplineGenerate(SuperResolution);
    
    id.loadSample(SuperResolution);
    PRINT(MSG_NORMAL, "Load sample!\n");
    is_color_image = false;////
    PRINT(MSG_NORMAL, "Is color image? %s\n", is_color_image ? "Yes" : "No");
    dm = new DataManager();
    dm->setWidth(id.width*SuperResolution);//sub-pixel
    dm->setHeight(id.height*SuperResolution);
    
    dm->initCUDA();
    dm->setData(id.getImage());
    clear_color = id.clear_color;
    dm->set_clear_color(id.clear_color);
}

void initShader() {
    string runifs[] = {"alpha",
                       "layer"
                      };
    sh_render = new GLSLProgram("glsl/render.vert", "glsl/render.frag");
    sh_render->compileAttachLink();
    sh_render->bindUniforms(ARRAY_SIZE(runifs), runifs);
}

void draw_image(DataManager* a_dm) {
    Texture2D* alpha = 0;
    glClear(GL_COLOR_BUFFER_BIT);
    
    alpha = a_dm->getAlphaMapOfLayer(impmap,5,SuperResolution);
            
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    /* SECOND PASS: Render using alpha map */
    sh_render->bind();
    glActiveTexture(GL_TEXTURE0);
    alpha->bind();
    //glViewport(0, 0, 1920, 1080);
    glUniform1f(sh_render->uniform("layer"), 1);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 1.0);
    glVertex2f(0, 0);
    glTexCoord2f(0.0, 0.0);
    glVertex2f(0, WHEIGHT);
    glTexCoord2f(1.0, 0.0);
    glVertex2f(WWIDTH, WHEIGHT);
    glTexCoord2f(1.0, 1.0);
    glVertex2f(WWIDTH, 0);
    glEnd();
    sh_render->unbind();
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    
}

void display(void) {
    dm->setClearColor();
    if (is_color_image) {
        dm_b->setClearColor();
        dm_g->setClearColor();
    }

    draw_image(dm);
    
    //glFlush();  // Finish rendering
    //glutSwapBuffers();

    if (SaveMore) //wang.
    {
        stringstream ss;
        ss<<"output"<<i<<".png";
        outFile = const_cast<char*>(ss.str().c_str());
    }
    else  outFile = const_cast<char*>("output.png");
    
    saveOutput();
    glutDestroyWindow(display1);
}

void reshape(int w, int h) {
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, WWIDTH, WHEIGHT, 0);
    glViewport(0, 0, WWIDTH, WHEIGHT);
    //glViewport(0, 0, 1920, 1080);
    glMatrixMode(GL_MODELVIEW);
    return;
}

void idle(void) {
    usleep(1000);
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
    case 's':
        outFile = const_cast<char*>("output.png");
        saveOutput();
        return;
        break;
    case 'q':
    case 27:
        exit(0);
        break;
    case 32:
        // TODO(maarten): interpolated image caching?
        interpolate = !interpolate;
        break;
    default:
        printf("Unsupported input: %c", key);
        fflush(stdout);
        break;
    }
    glutPostRedisplay();
}

int main(int argc, char *argv[]) {
    /* Read meta data and set variables accordingly */
    if (argc==4) //need to save lots of output images.
    {
        SaveMore = 1;
        i = argv[3];
    } /**/
    if (argc==3) 
    {
        SuperResolution = atoi (argv[2]);
    } 
    initDataManager(argv[1]);
   
    WHEIGHT = dm->getHeight();
    WWIDTH = dm->getWidth();
    PRINT(MSG_VERBOSE, "Image dimensions: %d x %d\n", WWIDTH, WHEIGHT);


    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(WWIDTH, WHEIGHT);
    
    display1 = glutCreateWindow(WINDOW_TITLE);

    // Initialize GLEW
    GLenum err = glewInit();
    if (GLEW_OK != err) {
        // Problem: glewInit failed, something is seriously wrong.
        PRINT(MSG_ERROR, "Error: %s\n", glewGetErrorString(err));
        exit(-1);
    }
    PRINT(MSG_NORMAL, "GLEW: Using GLEW %s\n", glewGetString(GLEW_VERSION));

    // Set OPENGL states 
    glEnable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // Set clear color one below the first layer that has points

    // Set texture parameters 
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);


    glutDisplayFunc(display);
     
   // glutKeyboardFunc(keyboard);
    glutReshapeFunc(reshape);
   // glutIdleFunc(idle);


    initShader();
    dm->initFBO();
    dm->setAlphaShader(SHADER);
    if (is_color_image) {
        dm_g->initFBO();
        dm_g->setAlphaShader(SHADER);
        dm_b->initFBO();
        dm_b->setAlphaShader(SHADER);
    }

    
    glutMainLoop();
    return 0;
}