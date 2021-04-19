#include "include/texture2D.hpp"
#include "include/vis.h"
#include "include/field.h"
#include "include/tinycolormap.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#define colorNum 1000
using namespace std;
typedef unsigned char BYTE;
	
Display*	Display::instance = 0;
float maxDistance, maxskel;
BYTE r,g,b;
int display1;
int DistinguishColor[colorNum][3];//wang

void generateColor()//wang
{
    srand (time(NULL));
    for (int i = 0; i<colorNum; i++)
    {
        DistinguishColor[i][0] = rand()%256;//generate a random num which lies in 0-255.
        DistinguishColor[i][1] = rand()%256;
        DistinguishColor[i][2] = rand()%256;
    }
}

Display::Display(int winSize_,int texsize, int argc,char** argv, FIELD<float>* splineField, FIELD<float>* boundColor, FIELD<float>* branchColor,FIELD<float>* CPmap, float hausdorff):
				imgSize(texsize),winSize(winSize_),scale(1),length(100),
				transX(0),transY(0),isLeftMouseActive(false),isRightMouseActive(false),oldMouseX(0),oldMouseY(0),
				show_what(0),tex_interp(true), splineF(splineField), bdColor(boundColor), bdColor2(branchColor),cpMap(CPmap)
{ 
    maxDistance = hausdorff;
	instance = this;    
    glutInitWindowSize(winSize, winSize); 
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA); 
    glutInit(&argc, argv); 
    display1 = glutCreateWindow("Boundary-rainbow colorspace"); 
    glutDisplayFunc(display);
	glutMouseFunc(mouse);
    glutMotionFunc(motion);
	glutKeyboardFunc(keyboard);
    generateColor();
    glGenTextures(1, &texture);
}

void Display::saveOutput() {
    char *outFile = const_cast<char*>("output.png");
    
    //int weight = 700; int height = instance->imgSize;
    int weight = 1000; int height = 900;
    unsigned char *sdata = (unsigned char *) malloc(weight * height * 4);
    glReadPixels(0, (1000-height), weight, height, GL_RGB, GL_UNSIGNED_BYTE, sdata);
    
    Texture2D *t = new Texture2D(weight, height, GL_RGB, GL_RGB, GL_UNSIGNED_BYTE);
    
    t->setData(sdata);
    t->saveAsPNG(outFile);
    free(sdata);
    delete t;
    /**/
}



void Display::display() 
{
    // Initialization
    glViewport(0, 0, (GLsizei) instance->winSize, (GLsizei) instance->winSize); 

    glClearColor(1.0, 1.0, 1.0, 1.0); 
    glClear(GL_COLOR_BUFFER_BIT); 

    // Setup projection matrix
    glMatrixMode(GL_PROJECTION); 
    glLoadIdentity(); 
    gluOrtho2D(0.0, 1.0, 0.0, 1.0); 

    // Setup modelview matrix
    glMatrixMode(GL_MODELVIEW); 
    glLoadIdentity(); 
    glScalef(instance->scale, instance->scale, 1.0);
    glTranslatef(instance->transX, instance->transY, 0.0);

    // Setting up lights
    glDisable(GL_LIGHTING); 
    glDisable(GL_DEPTH_TEST); 
    glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (instance->tex_interp)?GL_LINEAR:GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (instance->tex_interp)?GL_LINEAR:GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);


    // Draw the sample scene
    glBindTexture(GL_TEXTURE_2D, instance->texture); 


    glColor3f(1,1,1);
    glBegin(GL_QUADS);
    
    glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 1.0);
    glTexCoord2f(1.0, 0.0);  glVertex2f(1.0, 1.0);
    glTexCoord2f(1.0, 1.0);  glVertex2f(1.0, 0.0); 
    glTexCoord2f(0.0, 1.0);  glVertex2f(0.0, 0.0);	
    glEnd(); 
    
    saveOutput();//wang

    glDisable(GL_TEXTURE_2D); 
	glFinish();

    glutSwapBuffers();     

    glutDestroyWindow(display1);
}

void Display::float2rgb(float p)
{
    float v;
    if (!show_what)
        v= p/maxDistance;
    else
        v = p/maxskel;
        
    v = (v>0)?v:0; 
    v = (v<1)?v:1;
/*
//Viridis colormap    
    v = v*2.0/3.0; //cut the brightest third of the viridis colormap

    const tinycolormap::Color color = tinycolormap::GetColor(v, tinycolormap::ColormapType::Viridis);
    r = (BYTE)(255 * color.r());
    g = (BYTE)(255 * color.g());
    b = (BYTE)(255 * color.b());
*/
 //rainbow colormap
    const float dx=0.8f;

    v = (6-2*dx)*v+dx;
    r = (BYTE)(255 * max(0.0f,(3-(float)fabs(v-4)-(float)fabs(v-5))/2));
    g = (BYTE)(255 * max(0.0f,(4-(float)fabs(v-2)-(float)fabs(v-4))/2));
    b = (BYTE)(255 * max(0.0f,(3-(float)fabs(v-1)-(float)fabs(v-2))/2));

}

void Display::generateTexture()
{
    //bdColor->writePGM("bd.pgm");
	
    BYTE* tex = new BYTE[imgSize * imgSize * 3];			// Local buffer to store the GL texture
    
    for (int i = 0; i < imgSize; ++i)	//!!!!!!!!!!!!			// Generate visualization texture
        for (int j = 0; j < imgSize; ++j) 
		{
            int id = j * imgSize + i;															// Position of (i,j) in the various 1D arrays
 			
            if (!show_what) 
            {
                r=g=b=255;//white background
                /* */ 
                float* sf = splineF->data()+id;
                //cout<<"i: "<<i<<" j: "<<j<<" sf "<<*sf<<endl;
                if (*sf==128)
                    //r = g = b = 255; 
                    {r = 135; g = 206; b = 235; }//sky blue
                //else if (*sf==255)
                //    r = g = b = 0;
                
                float* p = bdColor->data()+id;
                if (*p)									// Boundary (site): mark as black.
                    { r = g = b = 0;}
               
                float* p2 = bdColor2->data()+id;
                int index = (int)(*p2);
                if (index){																
                     r = DistinguishColor[index][0];
                     g = DistinguishColor[index][1];
                     b = DistinguishColor[index][2];
                    }
                
                float* p3 = cpMap->data()+id;
                index = (int)(*p3);
                if (index){																
                     r = DistinguishColor[index][0];
                     g = DistinguishColor[index][1];
                     b = DistinguishColor[index][2];
                     //cout<<(int)r<<" "<<(int)g<<" "<<(int)b<<" ";
                    }/**/
                /*  
                if (*p)																				// Boundary (site): mark as red
                    float2rgb(*p);
                    //{ r = 255; g = b = 0;}
             
                 */  
            }
            
            tex[id * 3 + 0] = r;
            tex[id * 3 + 1] = g;
            tex[id * 3 + 2] = b; 
		}

    // Create the texture; the texture-id is already allocated
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imgSize, imgSize, 0, GL_RGB, GL_UNSIGNED_BYTE, tex); 

    delete[] tex;
}


void Display::mouse(int button, int state, int x, int y) 
{
    if (state == GLUT_UP)
        switch (button)
    {
        case GLUT_LEFT_BUTTON:
            instance->isLeftMouseActive = false;
            break;
        case GLUT_RIGHT_BUTTON:
            instance->isRightMouseActive = false; 
            break; 
    }

    if (state == GLUT_DOWN)
    {
        instance->oldMouseX = x;
        instance->oldMouseY = y;

        switch (button)
        {
        case GLUT_LEFT_BUTTON:
            instance->isLeftMouseActive = true; 
            break;
        case GLUT_RIGHT_BUTTON:
            instance->isRightMouseActive = true;
            break;
        }
    }		
}



void Display::motion(int x, int y) 
{
    if (instance->isLeftMouseActive) {
        instance->transX += 2.0 * double(x - instance->oldMouseX) / instance->scale / instance->imgSize; 
        instance->transY -= 2.0 * double(y - instance->oldMouseY) / instance->scale / instance->imgSize; 
        glutPostRedisplay(); 
    }
    else if (instance->isRightMouseActive) {
        instance->scale -= (y - instance->oldMouseY) * instance->scale / 400.0;
        glutPostRedisplay(); 
    } 

    instance->oldMouseX = x; instance->oldMouseY = y; 
}



void Display::keyboard(unsigned char k,int,int)
{
  int imgSize = instance->imgSize;
  int xm = instance->xm, ym = instance->ym, xM = instance->xM, yM = instance->yM;
  
  switch (k)
  {
    case '.':  instance->scale *= 0.9; break;
	case ',':  instance->scale *= 1.1; break;
	case ' ':  instance->show_what++; if (instance->show_what>1) instance->show_what=0; instance->generateTexture(); break;
	case 't':  instance->tex_interp = !instance->tex_interp; break;
	case 27:   exit(0);
  }
  
  glutPostRedisplay();
}





