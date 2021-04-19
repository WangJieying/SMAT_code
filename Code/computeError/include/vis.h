#pragma once


//	Display:		Very simple class to hold the visualization params (thus, not part of the skeletonization/FT engine per se)
//
//


#include <stdlib.h>
#include <GL/glut.h>
#include "include/field.h"


class	Display
{
public:


		Display(int winSize,int texsize, int argc, char** argv, FIELD<float>* splineField, FIELD<float>* boundColor, FIELD<float>* boundColor2,FIELD<float>* CPmap, float hausdorff
				);
void	generateTexture();
//void 	float2rgb(float p, BYTE r, BYTE g, BYTE b);
void 	float2rgb(float p);
//void DistinguishColor(int n);

int    imgSize;												//Effective size (pow 2) of texture images to be displayed

GLuint texture;												//Visualization parameters			
short* FT;
unsigned char* skel;
float* skelDT;
float* siteParam;
int    winSize; 
float  scale, transX, transY;
bool   isLeftMouseActive, isRightMouseActive; 
int    oldMouseX, oldMouseY;
int    show_what;
float  threshold;
short* endpoints;
int    nendpoints;
bool   tex_interp;
int    xm,ym,xM,yM;
const float length;
FIELD<float>* splineF;
FIELD<float>* bdColor, * bdColor2, *cpMap;
//FIELD<float>* skColor, * skColor2;


private:


static void mouse(int button, int state, int x, int y);
static void motion(int x, int y);
static void display();
static void saveOutput();
static void keyboard(unsigned char k,int,int);

static Display*	
			instance;
};






