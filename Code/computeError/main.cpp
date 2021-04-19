
#include "include/field.h"
#include "include/vis.h"
#include <cuda_runtime_api.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <time.h>

using namespace std;
//wang

#define INDEX(i,j) (i)+fboSize*(j)
typedef std::pair<int, int> coord2D_t;
typedef vector<coord2D_t> coord2D_list_t;
ofstream OutFile;

FIELD<float>* field, * field1, * deltafield, * origfield, * origfield1;
int width, height;
FIELD <float>* boundary, * boundary2, *origboundary_draw;
FIELD <float>* orig_boundary, * spline_boundary;
FIELD <float>* skel, * skel2, *skel_draw;
FIELD <float>* white_skel, * white_skel2, *white_skel_draw;
int fboSize;


int skelft2DSize(int nx,int ny)
{
	int dx = pow(2.0f,int(log(float(nx))/log(2.0f)+1));					//Find minimal pow of 2 which fits the input image
	int dy = pow(2.0f,int(log(float(ny))/log(2.0f)+1));	
	int fboSize	= max(dx,dy);											//DT/FT/skeleton image size (pow of 2, should be able to contain the input image)	
	
	return fboSize;
}


void convert(FIELD<float>* field)
{
    float *c = field->data();
    
    //float *end = field->data() + field->dimX*field->dimY;
    for (int y = 0; y < field->dimY(); ++y) 
        for (int x = 0; x < field->dimX(); ++x) 
        field->set(x, y, (255.0 - *c++));
    
    //field->writePGM("Orig.pgm");  
}

/**/
vector<coord2D_t> saveBoundary(FIELD<float>* field, int which)
{
    if (which == 1) orig_boundary = new FIELD<float> (fboSize, fboSize);
    if (which == 2)  spline_boundary = new FIELD<float> (fboSize, fboSize);

    vector<coord2D_t> Bound;

    for (int y = 0; y < field->dimY(); ++y) 
        for (int x = 0; x < field->dimX(); ++x) 
            if (!(*field)(x, y)) 
            {
                if(x==0||y==0||x==(field->dimX()-1)||y==(field->dimY()-1)) //lie in the boundary of the image.
                { 
                    Bound.push_back(coord2D_t(x,y));
                    if (which == 1) orig_boundary->set(x,y,1);
                    if (which == 2) spline_boundary->set(x,y,1);
                }
                else
                {
                    if((*field)(x-1, y)||(*field)(x+1, y)||(*field)(x, y-1)||(*field)(x, y+1))
                    {
                        Bound.push_back(coord2D_t(x,y));
                        if (which == 1) orig_boundary->set(x,y,1);
                        if (which == 2) spline_boundary->set(x,y,1);
                    }
                }
                 
            }
    return Bound;
}

vector<coord2D_t> saveSkel(FIELD<float>* field, bool first)
{
    if (first) white_skel = new FIELD<float> (fboSize, fboSize);
    else white_skel2 = new FIELD<float> (fboSize, fboSize);
    vector<coord2D_t> origBound;

    for (int y = 0; y < field->dimY(); ++y) 
        for (int x = 0; x < field->dimX(); ++x) 
            if ((*field)(x, y)) 
            {
                origBound.push_back(coord2D_t(x,y));
                if (first) white_skel->set(x,y,1);
                else white_skel2->set(x,y,1);
            }
    return origBound;
}

float countHausdorff(vector<coord2D_t> origBound, vector<coord2D_t> compBound, int whichone)
{
    cout<<"--"<<origBound.size()<<" "<<compBound.size()<<endl;
    switch (whichone)
    {
    case 1:
        boundary = new FIELD<float> (fboSize, fboSize); break; //orig-spline
    case 2:
        boundary2 = new FIELD<float> (fboSize, fboSize); break;//spline-orig
    case 3:
        skel = new FIELD<float> (fboSize, fboSize); break;
    case 4:
        skel2 = new FIELD<float> (fboSize, fboSize); break;
    }
    
        
    float maxLength = 0.0f;
    float sqrLength, minLength;
    float sum= 0;
    int count = 0;
    int x0,x1,y0,y1;
    vector<coord2D_t>::iterator it, itc;
    for(it=origBound.begin();it!=origBound.end();it++)
    {
        //cout<<(origBound.end()-it)<<"\t";
        x0 = (*it).first;
        y0 = (*it).second;
        minLength = 10000.0f;
        for(itc = compBound.begin();itc!=compBound.end();itc++)
        {
            x1 = (*itc).first;
            y1 = (*itc).second;
            sqrLength = (float)((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
            if (sqrLength < minLength) minLength = sqrLength;
            if (sqrLength == 0) { /**/minLength = 0.000001; break;}
        }

        switch (whichone)
        {
        case 1:
            boundary->set(x0,y0,sqrt(minLength)); break;//for the color bound.
        case 2:
            boundary2->set(x0,y0,sqrt(minLength)); break;//spline-orig
        case 3:
            skel->set(x0,y0,sqrt(minLength)); break;
        case 4:
            skel2->set(x0,y0,sqrt(minLength)); break;
        }

        //////support for the average H/////////
        sum += sqrt(minLength);
        count++;

        if (minLength > maxLength) maxLength = minLength;
        
        //if((int)(origBound.end()-it)<4) break;
    }
    cout<<sum<<" "<<count<<endl;
    OutFile<<"Average hausdorff: "<<sum/(float)count/sqrt(width*width+height*height)<<endl;
    return sqrt(maxLength);
}

void CountJaccard(FIELD<float>* f1, FIELD<float>* f2)
{
    int s1, s2, s12;//int is enough
    s1=s2=s12=0;
    float Jindex;
    for (int y = 0; y < f1->dimY(); ++y) 
        for (int x = 0; x < f1->dimX(); ++x) 
        {
            if (!(*f1)(x, y)) s1++;
            if (!(*f2)(x, y)) s2++;
            if ( (!(*f1)(x, y)) && (!(*f2)(x, y)) ) s12++;
        }
    cout<<s1<<" "<<s2<<" "<<s12<<endl;        
    Jindex = (float)s12 / (float)(s1+s2-s12);
    OutFile<<"Jaccard index "<<endl;
    OutFile<<Jindex<<endl;
}
/**/


int main(int argc,char **argv)
{
    
	if (argc<2) return 1;

	const char* filename = "../imShow/output.pgm";
    
    const char* origfile = argv[1];
   
    field = FIELD<float>::read(filename);
   
    origfield = FIELD<float>::read(origfile);
   
    if (!field||!origfield) {
        cout<<"Failed to read file."<<endl;
        exit(EXIT_FAILURE);
    }
    width = field->dimX();
    height = field->dimY();
    fboSize = skelft2DSize(width,height);//the width and height of the texture need to be the pow of 2.

    OutFile.open("../output.txt",ios_base::app);
    OutFile<<"================="<<endl;
    
////////////////calculate the boundary hausdorff////////////////
    convert(field); //color flip
    
    convert(origfield);
    vector<coord2D_t> origBound = saveBoundary(origfield, 1);
    vector<coord2D_t> compBound = saveBoundary(field, 2);

    float hausdorff = countHausdorff(origBound, compBound, 1);
    OutFile<<"Boundary hausdorff-(orig, spline): "<<endl;
    OutFile<<hausdorff<<endl;
    OutFile<<"Boundary hausdorff-(%): "<<hausdorff/sqrt(width*width+height*height)<<endl;
    OutFile<<"-----------------"<<endl;
    //float hausdorff2 = countHausdorff(compBound, origBound, 2);
   // OutFile<<"Boundary hausdorff-(spline, orig): "<<hausdorff2<<endl;
    //OutFile<<"================="<<endl;

    CountJaccard(field,origfield);
    OutFile<<"================="<<endl;
    
    OutFile<<"MS-SSIM: "<<endl;
    OutFile.close();


    field = FIELD<float>::read(filename);
   
    FIELD <float>* splinefield= new FIELD<float> (fboSize, fboSize);
     for (int y = 0; y < field->dimY(); ++y) 
        for (int x = 0; x < field->dimX(); ++x) 
        {
            if ((*field)(x, y)) splinefield->set(x,y,128);
            else splinefield->set(x,y,255);
        }
/**/
//////////read sample file-draw colorful branches.////
    FIELD <float>* branchColor= new FIELD<float> (fboSize, fboSize);
     
    string str;
    ifstream ifs1("BranchSample.txt"); 
    int x,y,v;
    while(ifs1)
    { 
        //count1++;
        ifs1 >> str;
        v = (int)atof(str.c_str());
        ifs1 >> str;
        x = (int)atof(str.c_str());
        ifs1 >> str;
        y = (unsigned int)atof(str.c_str());
        if(ifs1.fail())  break;

        branchColor->set(x,y,v);
        branchColor->set(x-1,y,v);
        branchColor->set(x+1,y,v);
        branchColor->set(x,y-1,v);
        branchColor->set(x,y+1,v);
 
    }
    /////////draw control points/////
    FIELD <float>* CPmap= new FIELD<float> (fboSize, fboSize);
    int branchN = 0;
    ifstream ifs2("ControlPoints.txt"); 
    //ifstream ifs2("1.txt"); 
    ifs2 >> str;
    ifs2 >> str;
    
    int CPnum;
    while(ifs2)
    { 
        branchN++;
        ifs2 >> str;
        CPnum = (int)atof(str.c_str());
        ifs2 >> str;
        //degree = (int)atof(str.c_str());
        ifs2 >> str;
        //numSamples = (unsigned int)atof(str.c_str());
    
        if(ifs2.fail())  break;
        
        for (int i = 0; i< CPnum; ++i)
        {
            
            ifs2 >> str;
            x = (int)atof(str.c_str());
            ifs2 >> str;
            y = (int)atof(str.c_str());
            ifs2 >> str;
            
            //draw thicker
            for(int k =-3;k<4;k++)
                CPmap->set(x+k,y,branchN);
            for(int k =-2;k<3;k++)
                CPmap->set(x+k,y-1,branchN);
            for(int k =-2;k<3;k++)
                CPmap->set(x+k,y+1,branchN);
            for(int k =-1;k<2;k++)
                CPmap->set(x+k,y-2,branchN);
            for(int k =-1;k<2;k++)
                CPmap->set(x+k,y+2,branchN);
            CPmap->set(x,y-3,branchN);
            CPmap->set(x,y+3,branchN); 
        }
    }
    ///end
    Display* dpy = new Display(1000,fboSize,argc,argv, splinefield, boundary,branchColor,CPmap,hausdorff);
																		//Initialize visualization engine
    dpy->generateTexture();												//Show results
    glutMainLoop(); 	

    //deinitialization(siteParam, outputFT, outputSkeleton, outputTopo, skelDT);
	delete dpy;
    return 0;
}



