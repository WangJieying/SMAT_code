/* skeleton.cpp */

#include <math.h>
#include "fileio/fileio.hpp"
#include <cuda_runtime_api.h>
#include "CUDASkel2D/include/genrl.h"
#include "CUDASkel2D/include/field.h"
#include "CUDASkel2D/include/skelft.h"
#include "ImageWriter.hpp"
#include "skeleton_cuda.hpp"
#include "include/connected.hpp"
#include "include/messages.h"


float SKELETON_SALIENCY_THRESHOLD;
float IMP_THRESHOLD;
float SKELETON_ISLAND_THRESHOLD;
float        SKELETON_DT_THRESHOLD;

#define INDEX(i,j) (i)+fboSize*(j)
 
float* siteParam;
// The following are CUDA buffers for several images.
//unsigned char* outputSkeleton;
float* outputSkeleton;
short* outputFT;
short* skelFT;
//float* outputDT;
bool* foreground_mask;
int xm = 0, ym = 0, xM, yM, fboSize;
int dimX, dimY;
double length;

void allocateCudaMem(int size) {
    skelft2DInitialization(size);
    cudaMallocHost((void**)&outputFT, size * size * 2 * sizeof(short));
    cudaMallocHost((void**)&foreground_mask, size * size * sizeof(bool));
    cudaMallocHost((void**)&outputSkeleton, size * size * sizeof(float));
    cudaMallocHost((void**)&siteParam, size * size * sizeof(float));
    cudaMallocHost((void**)&skelFT, size * size * 2 * sizeof(short));
    //cudaMallocHost((void**)&outputDT,size*size*2*sizeof(float));//wang
}

void deallocateCudaMem() {
    skelft2DDeinitialization();
    cudaFreeHost(outputFT);
    cudaFreeHost(foreground_mask);
    cudaFreeHost(outputSkeleton);
    cudaFreeHost(siteParam);
    cudaFreeHost(skelFT);
    //cudaFreeHost(outputDT);
}


int initialize_skeletonization(FIELD<float>* im) {
    xM = im->dimX();
    yM = im->dimY();
    dimX = im->dimX();
    dimY = im->dimY();
    fboSize = skelft2DSize(xM, yM);//   Get size of the image that CUDA will actually use to process our nx x ny image

    allocateCudaMem(fboSize);
    return fboSize;
}

FIELD<float>* skel_to_field(FIELD<float>* im) {
    FIELD<float>* f = new FIELD<float>(dimX, dimY);
    for (int i = 0; i < dimX; ++i) {
        for (int j = 0; j < dimY; ++j) {
            //bool is_skel_point = outputSkeleton[INDEX(i, j)];
            //f->set(i, j, is_skel_point ? 255 : 0);
            f -> set(i,j,outputSkeleton[INDEX(i, j)]);
            //cout<<outputSkeleton[INDEX(i, j)]<<" ";
            //if(im->value(i,j)==0 && (im->value(i-1,j) || im->value(i+1,j) || im->value(i,j-1) || im->value(i,j+1))) f -> set(i,j,0);
        }
    }
    return f;
}
/*
void getDT(FIELD<float>* f)  //why does "FIELD<float>* getDT()" cannot change input value?
{
    for (int i = 0; i < xM; ++i) {
        for (int j = 0; j < yM; ++j) {
            f->set(i, j, outputDT[INDEX(i,j)]);
            //cout<<outputDT[INDEX(i,j)]<<" ";
        }
    }
}
*/

// TODO(maarten): dit moet beter kunnen, i.e. in een keer de DT uit cuda halen
void dt_to_field(FIELD<float>* f) {
    for (int i = 0; i < xM; ++i) {
        for (int j = 0; j < yM; ++j) {
            int id = INDEX(i, j);
            int ox = outputFT[2 * id];
            int oy = outputFT[2 * id + 1];
            float val = sqrt((i - ox) * (i - ox) + (j - oy) * (j - oy));
            f->set(i, j, val);
        }
    }
}

FIELD<float>* skelft_to_field() {
 
    FIELD<float>* f = new FIELD<float>(xM, yM);
    
    for (int i = 0; i < xM; ++i) {
        for (int j = 0; j < yM; ++j) {
            int id = INDEX(i, j);
            int ox = skelFT[2 * id];
           // PRINT(MSG_NORMAL, "f2: %d\n", ox);
            int oy = skelFT[2 * id + 1];
            double val = 255 * (siteParam[oy * fboSize + ox] / length);
            f->set(i, j, val);
        }
    }
    return f;
}

void analyze_cca(FIELD<float>* skel) {
   // skel -> writePGM("before.pgm");
    float* ssd = skel->data();
    
    ConnectedComponents *CC = new ConnectedComponents(255);
    int                 *ccaOut = new int[skel->dimX() * skel->dimY()];
    int                 hLabel;
    unsigned int        *hist;
    int                 i,
                        nPix = skel->dimX() * skel->dimY();

    /* Perform connected component analysis */
    hLabel = CC->connected(ssd, ccaOut, skel->dimX(), skel->dimY(), std::equal_to<float>(), true);
    hist = static_cast<unsigned int*>(calloc(hLabel + 1, sizeof(unsigned int)));    
    if (!hist) {
        PRINT(MSG_ERROR, "Error: Could not allocate histogram for skeleton connected components\n");
        exit(-1);
    }
    //for (i = 0; i < hLabel; i++) { hist[i] = 0; /*if (hist[i] != 0) {/*cout<<hist[i]<<endl; counNonZero =+ 1;}*/ }
 
    for (i = 0; i < nPix; i++) { hist[ccaOut[i]]++; }
    
    /* Remove small islands */
    for (i = 0; i < nPix; i++) {
         ssd[i] = (hist[ccaOut[i]] >= SKELETON_ISLAND_THRESHOLD/100*sqrt(skel->dimX()*skel->dimY())) ? ssd[i] : 0;
      // ssd[i] = (hist[ccaOut[i]] >= SKELETON_ISLAND_THRESHOLD && hist[ccaOut[i]] <= SKELETON_ISLAND_THRESHOLD * 2) ? ssd[i] : 0;
        
       // ssd[i] = (hist[ccaOut[i]] >= SKELETON_ISLAND_THRESHOLD) ? ssd[i] : 0;
    }

    free(hist);
    delete CC;
    delete [] ccaOut;

    skel->setData(ssd);
   // skel->writePGM("removed.pgm");
}

FIELD<float>* computeSkeleton(FIELD<float> *input) {
   
    memset(siteParam, 0, fboSize * fboSize * sizeof(float));
   
    memset(foreground_mask, false, fboSize * fboSize * sizeof(bool));
    
    unsigned char* image = new unsigned char[fboSize*fboSize];


    int nx = input->dimX();
    int ny = input->dimY();
    xm = ym = nx; xM = yM = 0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (!(*input)(i, j)) {
                foreground_mask[INDEX(i, j)] = true; //original
                //if(!(*input)(i-1, j) && !(*input)(i+1, j) && !(*input)(i, j-1) && !(*input)(i, j+1)) foreground_mask[INDEX(i, j)] = true; //wang. Remove the boundary
                siteParam[INDEX(i, j)] = 1;
                image[j*fboSize+i] = 255;
                xm = min(xm, i); ym = min(ym, j);
                xM = max(xM, i); yM = max(yM, j);
            }
        }
    }
    
    xM = nx; yM = ny;//orginal method, which cause artifacts.
    
   /* 
    xM = (xm==nx)? (nx - 1): nx;//must run before xm=...
    yM = (xm==nx)? (ny - 1): ny;
    xm = (xm==0)? 0 : (xm-1); // In this way, we further expand the boundary, which avoid producing artifacts.
    ym = (ym==0)? 0 : (ym-1);
   */
    skelft2DFT(0, siteParam, 0, 0, fboSize, fboSize, fboSize);
    skelft2DDT(outputFT, 0, xm, ym, xM, yM); 
   // length = skelft2DMakeBoundary((unsigned char*)outputFT, xm, ym, xM, yM, siteParam, fboSize, 0, false);
  
    length = skelft2DMakeBoundary((unsigned char*)outputFT, xm, ym, xM, yM, siteParam, fboSize, 1, true);
   // length = skelft2DMakeBoundary((unsigned char*)image,xm,ym,xM,yM,siteParam,fboSize,1,true);
	/* //The following works the same as the one in image, so I keep that one.
    FIELD<float>* outdt = new FIELD<float> (*input);
    int count = 0;
    for (int i = 0; i < xM; ++i) {
        for (int j = 0; j < yM; ++j) {
            
            if (siteParam[INDEX(i, j)]>0)
                { outdt->set(i, j, 255.0); count++;}
            else
                outdt->set(i, j, 0.0);             
                
        }
    }
    cout<<"contour count: "<<count<<endl;
    
    //outdt -> writePGM("boundary.pgm");
    */
    if (!length) return NULL;
    skelft2DFillHoles((unsigned char*)outputFT, xm + 1, ym + 1, 1);
    skelft2DFT(outputFT, siteParam, xm, ym, xM, yM, fboSize);
    //cout<<foreground_mask[INDEX(172, 26)]<<"---"<<foreground_mask[INDEX(172, 27)]<<"---"<<foreground_mask[INDEX(171, 27)]<<endl;
    skelft2DSkeleton(outputSkeleton, foreground_mask, length, SKELETON_SALIENCY_THRESHOLD,IMP_THRESHOLD, xm, ym, xM, yM);
    
    skel2DSkeletonFT(skelFT, xm, ym, xM, yM);
   // skelft2DDT(outputDT,xm,ym,xM,yM);//wang
   // FIELD<float>* dt = getDT();
    //getDT(input);
    auto imp = skel_to_field(input);
    dt_to_field(input);
    //analyze_cca(skel);
    return imp;
}

short* get_current_skel_ft() {
    size_t sz = fboSize * fboSize * 2;
    short* skelft = (short*)(malloc(sz * sizeof(short)));
    if (!skelft) {
        PRINT(MSG_ERROR, "Error: could not allocate array for skeleton FT\n");
        exit(-1);
    }
    memcpy(skelft, skelFT, sz);
    return skelft;
}
