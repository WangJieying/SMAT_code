#include "include/Image.hpp"
#include <sys/stat.h>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <set>
#include "include/ImageWriter.hpp"
#include "include/skeleton_cuda.hpp"
#include "include/connected.hpp"
#include "include/messages.h"
#include <parallel/algorithm>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include "SplineGenerate/BSplineCurveFitter/BSplineCurveFitterWindow3.h"
#include <time.h> 


string OUTPUT_FILE;
bool BUNDLE, OVERLAP_PRUNE, MSIL_SELECTION;
float ALPHA;
float OVERLAP_PRUNE_EPSILON;
int EPSILON;
int MIN_PATH_LENGTH;
int MIN_SUM_RADIUS;
int MIN_OBJECT_SIZE;
//ofstream OutFile;
ofstream OutFile2;
bool firsttrace = true;
float ByteCount = 4;//width+height
FIELD<float>* skelImp;//for importance map
//FIELD<float>* segment;
vector<vector<Vector3<float>>> BranchSet;
float diagonal;
float hausdorff;
vector<int *> connection;
int SplitNum;
//int repeatNum = 0;

/*************** CONSTRUCTORS ***************/
Image::Image(FIELD<float> *in, float islandThresh, float importanceThresh) {
    PRINT(MSG_NORMAL, "Creating Image Object...\n");
    this->layerThreshold = importanceThresh;
    this->islandThreshold = islandThresh;
    this->importance = NULL;
    this->im = in;
    this->nPix = in->dimX() * in->dimY();
    std::set<int> levels(in->data(), in->data() + this->nPix);
    this->graylevels = reinterpret_cast<int*>(malloc(levels.size() * sizeof(int)));
    std::set<int>::iterator it = levels.begin();
    for (uint i = 0; i < levels.size(); ++i) {
        this->graylevels[i] = *it;
        std::advance(it, 1);
    }
    this->numLayers = levels.size();
    PRINT(MSG_NORMAL, "Done!\n");
}

Image::Image(FIELD<float> *in) {
    Image(in, 0, 0);
}

Image::~Image() {
    deallocateCudaMem();
    free(importance);
    free(graylevels);
}

/*************** FUNCTIONS **************/

/* rmObject -Remove current object in a 3x3 kernel, used for removeDoubleSkel: */
void rmObject(int *k, int x, int y) {
    if (x < 0 || x > 2 || y < 0 || y > 2 || k[y * 3 + x] == 0) { return; }
    k[y * 3 + x] = 0;
    rmObject(k, x + 1, y + 1);
    rmObject(k, x + 1, y);
    rmObject(k, x + 1, y - 1);
    rmObject(k, x, y + 1);
    rmObject(k, x, y - 1);
    rmObject(k, x - 1, y + 1);
    rmObject(k, x - 1, y);
    rmObject(k, x - 1, y - 1);
}
/* numObjects - Count the number of objects in a 3x3 kernel, used for removeDoubleSkel: */
int numObjects(int *k) {
    int c = 0;
    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; ++y) {
            if (k[y * 3 + x]) { rmObject(k, x, y); c++; }
        }
    }
    return c;
}
/* End count code */

/**
* removeDoubleSkel
* @param FIELD<float> * layer -- the layer of which the skeleton should be reduced
* @return new FIELD<float> *. Copy of 'layer', where all redundant skeleton-points are removed (i.e. rows of width 2.)
*/
void removeDoubleSkel(FIELD<float> *layer) {
    int *k = (int *)calloc(9, sizeof(int));
    for (int y = 0; y < layer->dimY(); ++y) {
        for (int x = 0; x < layer->dimX(); ++x) {
            if (layer->value(x, y)) {
                k[0] = layer->value(x - 1, y - 1);
                k[1] = layer->value(x - 1, y);
                k[2] = layer->value(x - 1, y + 1);
                k[3] = layer->value(x  , y - 1);
                k[4] = 0;
                k[5] = layer->value(x  , y + 1);
                k[6] = layer->value(x + 1, y - 1);
                k[7] = layer->value(x + 1, y);
                k[8] = layer->value(x + 1, y + 1);
                if (k[0] + k[1] + k[2] + k[3] + k[4] + k[5] + k[6] + k[7] + k[8] > 256) {
                    int b = numObjects(k);
                    if (b < 2) {layer->set(x, y, 0);}
                }
            }
        }
    }
    free(k);
    //return layer;
}

void NewRemoveDoubleSkel()
{
    //minimalSkel = new FIELD<float>(skelImp->dimX(), skelImp->dimY());
    int *k = (int *)calloc(9, sizeof(int));
    int thres = 8;
    for (int y = 0; y < skelImp->dimY(); ++y) {
        for (int x = 0; x < skelImp->dimX(); ++x) {
            if (skelImp->value(x, y)) {
                k[0] = skelImp->value(x - 1, y - 1);
                k[1] = skelImp->value(x - 1, y);
                k[2] = skelImp->value(x - 1, y + 1);
                k[3] = skelImp->value(x  , y - 1);
                k[4] = skelImp->value(x  , y);
                k[5] = skelImp->value(x  , y + 1);
                k[6] = skelImp->value(x + 1, y - 1);
                k[7] = skelImp->value(x + 1, y);
                k[8] = skelImp->value(x + 1, y + 1);
               
                if (k[0]) {// if (lower neighbor < 2 * lower left neighbor && lower neighbor < 2 * central point), then remove the lower neighbor; 
                    if(k[3] && k[3] < 2*k[0] && k[3] < 2*k[4]) skelImp->set(x, y-1, 0);
                    if(k[1] && k[1] < 2*k[0] && k[1] < 2*k[4]) skelImp->set(x-1, y, 0);
                }
                if(k[2]) {
                    if(k[1] && k[1] < 2*k[2] && k[1] < 2*k[4]) skelImp->set(x-1, y, 0);
                    if(k[5] && k[5] < 2*k[2] && k[5] < 2*k[4]) skelImp->set(x, y+1, 0);
                }
                if(k[6]) {
                    if(k[3] && k[3] < 2*k[6] && k[3] < 2*k[4]) skelImp->set(x, y-1, 0);
                    if(k[7] && k[7] < 2*k[6] && k[7] < 2*k[4]) skelImp->set(x+1, y, 0);
                }
                if(k[8]) {
                    if(k[7] && k[7] < 2*k[8] && k[7] < 2*k[4]) skelImp->set(x+1, y, 0);
                    if(k[5] && k[5] < 2*k[8] && k[5] < 2*k[4]) skelImp->set(x, y+1, 0);
                }
                /*
                if(k[5])
                    if(k[5] < thres) skelImp->set(x, y+1, 0);
                    else if(k[8] > thres1 && k[4] > thres1) skelImp->set(x, y+1, 0);
                */
            }
        }
    }
    free(k);
    //return layer;

}


coord2D_list_t * neighbours(int x, int y, FIELD<float> *skel) {
    coord2D_list_t *neigh = new coord2D_list_t();
    int n[8] = {1, 1, 1, 1, 1, 1, 1, 1};

    /* Check if we are hitting a boundary on the image */
    if (x <= 0 )             {        n[0] = 0;        n[3] = 0;        n[5] = 0;    }
    if (x >= skel->dimX() - 1) {        n[2] = 0;        n[4] = 0;        n[7] = 0;    }
    if (y <= 0)              {        n[0] = 0;        n[1] = 0;        n[2] = 0;    }
    if (y >= skel->dimY() - 1) {        n[5] = 0;        n[6] = 0;        n[7] = 0;    }

    /* For all valid coordinates in the 3x3 region: check for neighbours*/
    if ((n[0] != 0) && (skel->value(x - 1, y - 1) > 0)) { neigh->push_back(coord2D_t(x - 1, y - 1)); }
    if ((n[1] != 0) && (skel->value(x    , y - 1) > 0)) { neigh->push_back(coord2D_t(x    , y - 1)); }
    if ((n[2] != 0) && (skel->value(x + 1, y - 1) > 0)) { neigh->push_back(coord2D_t(x + 1, y - 1)); }
    if ((n[3] != 0) && (skel->value(x - 1, y    ) > 0)) { neigh->push_back(coord2D_t(x - 1, y    )); }
    if ((n[4] != 0) && (skel->value(x + 1, y    ) > 0)) { neigh->push_back(coord2D_t(x + 1, y    )); }
    if ((n[5] != 0) && (skel->value(x - 1, y + 1) > 0)) { neigh->push_back(coord2D_t(x - 1, y + 1)); }
    if ((n[6] != 0) && (skel->value(x    , y + 1) > 0)) { neigh->push_back(coord2D_t(x    , y + 1)); }
    if ((n[7] != 0) && (skel->value(x + 1, y + 1) > 0)) { neigh->push_back(coord2D_t(x + 1, y + 1)); }

    return neigh;
}
coord2D_list_t * is8not4(int x, int y, FIELD<float> *skel) {
    coord2D_list_t *neigh = new coord2D_list_t();
    int n[4] = {1, 1, 1, 1};

    /* Check if we are hitting a boundary on the image */
    if (x <= 0 )             {        n[0] = 0;        n[2] = 0;  }
    if (x >= skel->dimX() - 1) {        n[1] = 0;        n[3] = 0;  }
    if (y <= 0)              {        n[0] = 0;        n[1] = 0; }
    if (y >= skel->dimY() - 1) {        n[2] = 0;        n[3] = 0;  }

    /* For all valid coordinates in the 3x3 region: check for neighbours*/
    if ((n[0] != 0) && (skel->value(x - 1, y - 1) > 0)) { neigh->push_back(coord2D_t(x - 1, y - 1)); }
    if ((n[1] != 0) && (skel->value(x + 1, y - 1) > 0)) { neigh->push_back(coord2D_t(x + 1, y - 1)); }
    if ((n[2] != 0) && (skel->value(x - 1, y + 1) > 0)) { neigh->push_back(coord2D_t(x - 1, y + 1)); }
    if ((n[3] != 0) && (skel->value(x + 1, y + 1) > 0)) { neigh->push_back(coord2D_t(x + 1, y + 1)); }

    return neigh;
}


void tracePath(int x, int y, FIELD<float> *skel, FIELD<float> *dt, int &seq) {
    coord2D_t n, list1;
    coord2D_list_t *neigh, *neigh_;
    coord2D_list_t neighList; //vector<pair<int, int>>
    coord2D_list_t *SetZero = new coord2D_list_t();

    if (firsttrace)
    {
        firsttrace = false;

        FIELD<float>* FindEndSkel = new FIELD<float> (*skel);
        FindEndSkel->set(x, y, 0);
        neigh = neighbours(x, y, FindEndSkel);

        if (neigh->size()<2) ; //it's an end point.
        else
        {
            for (auto i = neigh->begin(); i!=neigh->end();i++)
                neighList.push_back(*i); 

            while(!neighList.empty()) 
            {
                list1 = neighList.front();
                //cout<<"list1.first: "<<list1.first<<" list1.second: "<<list1.second<<endl;
                neighList.erase(neighList.begin());//erase the first one.
                FindEndSkel -> set(list1.first, list1.second, 0);

                coord2D_list_t * neigh1 = neighbours(list1.first, list1.second, FindEndSkel);
                if (neigh1->size()==0) //end
                    {x = list1.first; y = list1.second; break;}
                else
                {
                    for (auto i = neigh1->begin(); i!=neigh1->end();i++)
                        neighList.push_back(*i);
                }
            }    
        }
    }
    
		//OutFile<<seq<<" " << x << " " << y << " " <<dt->value(x, y)<< " " <<skelImp->value(x,y) <<endl;//for 4D
    //OutFile<<seq<<" " << x << " " << y << " " <<dt->value(x, y)<<endl;//for 3D
	skel->set(x, y, 0);
    //segment->set(x,y,255);
	neigh = neighbours(x, y, skel);

	if (neigh->size() < 1) ;//==0; at end.
	else if (neigh->size() < 2)//==1; same branch.
	{
		n = *(neigh->begin());
		tracePath(n.first, n.second, skel, dt, seq);//same seq;
	}
	else if (neigh->size() > 1) {//produce branches.
		//seqT=seq;
		while (neigh->size() > 0) {
			seq++;
			
			n = *(neigh->begin());
            int first = n.first;
            int second = n.second;
            neigh->erase(neigh->begin());

            for (auto i = neigh->begin(); i!=neigh->end();i++){
                skel->set((*i).first,(*i).second, 0);//first set to 0
                if(abs(first-(*i).first) + abs(second-(*i).second)==1)//is nerghbor
                {
                    neigh_ = neighbours((*i).first, (*i).second, skel);
                    for (auto i_ = neigh_->begin();i_!=neigh_->end();i_++)
                    {
                        if ((*i_).first==first && (*i_).second==second) continue;
                        else
                        {
                            if(abs(skelImp->value((*i_).first,(*i_).second)-skelImp->value((*i).first,(*i).second))<5)
                             {
                                skel->set((*i_).first,(*i_).second,0);
                                SetZero->push_back(coord2D_t((*i_).first, (*i_).second));
                            }
                        }
                    }
                }
            }
                
			tracePath(first, second, skel, dt, seq);

            for (auto i = neigh->begin(); i!=neigh->end();i++)
                skel->set((*i).first,(*i).second, 1);//then set back to 1.

            for (auto i = SetZero->begin(); i!=SetZero->end();i++)
                skel->set((*i).first,(*i).second, 1);
            SetZero->clear();
			//delete neigh;
			//neigh = neighbours(x, y, skel);
    	}
	}
	 delete neigh;		
    
}

//Add a vector to store the connection.
void tracePath(int x, int y, FIELD<float> *skel, FIELD<float> *dt, vector<Vector3<float>> Branch, int &seq){
    coord2D_t list1; coord2D_list_t neighList;
    int *EachConnect = (int *)calloc(4, sizeof(int));
    int index = 0;
    Vector3<float> CurrentPx;
    coord2D_t n;
    coord2D_list_t *neigh, *neigh_;
    coord2D_list_t *SetZero = new coord2D_list_t();//vector<pair<int, int>>


    if (firsttrace)
    {
        firsttrace = false;

        FIELD<float>* FindEndSkel = new FIELD<float> (*skel);
        FindEndSkel->set(x, y, 0);
        neigh = neighbours(x, y, FindEndSkel);

        if (neigh->size()<2) ; //it's an end point.
        else
        {
            for (auto i = neigh->begin(); i!=neigh->end();i++)
                neighList.push_back(*i); 

            while(!neighList.empty()) 
            {
                list1 = neighList.front();
                //cout<<"list1.first: "<<list1.first<<" list1.second: "<<list1.second<<endl;
                neighList.erase(neighList.begin());//erase the first one.
                FindEndSkel -> set(list1.first, list1.second, 0);

                coord2D_list_t * neigh1 = neighbours(list1.first, list1.second, FindEndSkel);
                if (neigh1->size()==0) //end
                    {x = list1.first; y = list1.second; break;}
                else
                {
                    for (auto i = neigh1->begin(); i!=neigh1->end();i++)
                        neighList.push_back(*i);
                }
            }    
        }
    }



    CurrentPx[0] = (float)x/diagonal;
    CurrentPx[1] = (float)y/diagonal;
    CurrentPx[2] = (float)dt->value(x, y)/diagonal;

	//OutFile<<seq<<" " << x << " " << y << " " <<dt->value(x, y)<< " " <<skelImp->value(x,y) <<endl;//for 4D
   
    Branch.push_back(CurrentPx);
    
	skel->set(x, y, 0);
    //segment->set(x,y,255);
	neigh = neighbours(x, y, skel);

	if (neigh->size() < 1) BranchSet.push_back(Branch);//==0; at end.
	else if (neigh->size() < 2)//==1; same branch.
	{
		n = *(neigh->begin());
		tracePath(n.first, n.second, skel, dt, Branch, seq);//same seq;
	}
	else if (neigh->size() > 1) {//produce branches.
    
        BranchSet.push_back(Branch);
        EachConnect[index] = seq;
        //cout<<"seq: "<<seq<<endl;
        //repeatNum += neigh->size();
		while (neigh->size() > 0) {
			seq++;
            index++;
			n = *(neigh->begin());
            int first = n.first;
            int second = n.second;
            neigh->erase(neigh->begin());

            for (auto i = neigh->begin(); i!=neigh->end();i++){//See the explainaton in my note
                skel->set((*i).first,(*i).second, 0);//first set to 0  
                //open in important version //shut in saliency version
               /* if(abs(first-(*i).first) + abs(second-(*i).second)==1)//is nerghbor, then only set the point i to 0 is not enough, still need to check the neighbor of this point.
                {
                    neigh_ = neighbours((*i).first, (*i).second, skel);
                    for (auto i_ = neigh_->begin();i_!=neigh_->end();i_++)
                    {
                        if ((*i_).first==first && (*i_).second==second) continue;
                        else
                        {
                            if(abs(skelImp->value((*i_).first,(*i_).second)-skelImp->value((*i).first,(*i).second))<5)
                            {
                                skel->set((*i_).first,(*i_).second,0);
                                SetZero->push_back(coord2D_t((*i_).first, (*i_).second));
                            }
                        }
                    }
                
                } */
            }
            EachConnect[index] = seq;
            vector<Vector3<float>> NewBranch;
            NewBranch.push_back(CurrentPx);
			tracePath(first, second, skel, dt, NewBranch,seq);

            for (auto i = neigh->begin(); i!=neigh->end();i++)
                skel->set((*i).first,(*i).second, 1);//then set back to 1.

            for (auto i = SetZero->begin(); i!=SetZero->end();i++)
                skel->set((*i).first,(*i).second, 1);
            SetZero->clear();
			//delete neigh;
			//neigh = neighbours(x, y, skel);
    	}
        connection.push_back(EachConnect);
        //cout<<"here "<<connection.size()<<endl;
	}
	 delete neigh;	
}

skel_tree_t *tracePath(int x, int y, FIELD<float> *skel, FIELD<float> *dt) {
    coord2D_t n;
    coord2D_list_t *neigh;
    skel_tree_t *path;
    if (skel->value(x, y) == 0) {
        PRINT(MSG_ERROR, "Reached invalid point.\n");
        return NULL;
    }

    // Create new node and add to root
    path = new SkelNode<coord3D_t>(coord3D_t(x, y, dt->value(x, y)));
    skel->set(x, y, 0);

    neigh = neighbours(x, y, skel);
    // Add children
    while (neigh->size() > 0) {
        n = *(neigh->begin());
        path->addChild(tracePath(n.first, n.second, skel, dt));
        delete neigh;
        neigh = neighbours(x, y, skel);
    }

    delete neigh;
    return path;
}

skel_tree_t* traceLayer(FIELD<float>* skel, FIELD<float>* dt) {
    skel_tree_t* root;
    coord3D_t rootCoord = coord3D_t(-1, -1, -1);
    root = new skel_tree_t( rootCoord );
    for (int y = 0; y < skel->dimY(); ++y) {
        for (int x = 0; x < skel->dimX(); ++x) {
            if (skel->value(x, y) > 0) {
                root->addChild(tracePath(x, y, skel, dt));
            }
        }
    }
    return root;
}


void removeRtZero(FIELD<float>* skel, FIELD<float>* dt)
{
    for (int y = 0; y < skel->dimY(); ++y) {
        for (int x = 0; x < skel->dimX(); ++x) {
            if (skel->value(x, y) > 0) {
                if (dt->value(x, y)==0)
                skel->value(x, y)=0;
            }
        }
    }

}
//support for encodeforcontour() function.
void traceContour(int x, int y, FIELD<float> *skel) {
    coord2D_t n;
    coord2D_list_t *neigh;
     
    skel->set(x, y, 0);

    neigh = neighbours(x, y, skel);
    if (neigh->size()==0) //reach endpoint.
    {
        //cout<<"endpoint: "<<x <<" "<<y<<" "<<ByteCount<<endl;
        if (abs((round)(ByteCount) - ByteCount) > 0) // *.5 - means it needed to be make up.
            ByteCount += 0.5; // add 1111 to make it up.
        ByteCount += 1; //add 11111111 as the END mark.
    }
    else if (neigh->size() == 1)
        {
            n = *(neigh->begin());
            ByteCount += 0.5; //only needs 0.5 bytes for one point.
            //cout<<"1 "<<n.first <<" "<<n.second<<" "<<ByteCount<<"---";
            traceContour(n.first, n.second, skel);
        }
        else
        {
            int mind = 100, minx, miny; 
            for (auto i = neigh->begin(); i!=neigh->end();i++)
            {
                int d = abs(x - (*i).first) + abs(y - (*i).second);
                if (d < mind) {mind = d; minx = (*i).first; miny = (*i).second;}
            }

            ByteCount += 0.5; //only needs 0.5 bytes for one point.
            //cout<<"2 "<<minx <<" "<<miny<<" "<<ByteCount<<"---";
            traceContour(minx, miny, skel);   
        }
        
    
    delete neigh; 
}

/* 
* Calculate the space for storing the contour of a binary image (in two way).
*/
void Image::encodeforcontour()
{
    // Find the contour.
    OutFile2.open("../output.txt",ios_base::app);
   
    FIELD<float>* imDt = im->dupe();
    FIELD<float>* contour = new FIELD<float>(imDt->dimX(), imDt->dimY());
    int count = 0;
    for (int y = 0; y < imDt->dimY(); ++y) 
        for (int x = 0; x < imDt->dimX(); ++x) 
            if (imDt->value(x, y) == 0) 
			{
                if ( imDt->value(x-1, y) || imDt->value(x+1, y) || imDt->value(x, y-1) || imDt->value(x, y+1))
                    {contour -> set(x,y,1);  count++;}
                else if ( x==0 || y==0 || x==imDt->dimX()-1 || y== imDt->dimY()-1)
                        {contour -> set(x,y,1);  count++;}
            }
    cout<<"count "<<count<<endl;
    OutFile2<<"Storing contour pixels needs: "<<4 * count + 4<<" Bytes."<<endl;//the last 4 is used for width and height.
    
    //contour->writePGM("contour.pgm");
    removeDoubleSkel(contour);
    //cout<<"count "<<count<<endl;
    for (int y = 0; y < imDt->dimY(); ++y) 
        for (int x = 0; x < imDt->dimX(); ++x) 
            if (contour->value(x, y) > 0) 
			{
                ByteCount += 4; //needs store absolute position
                //cout<<"absolute: "<<x <<" "<<y<<" "<<ByteCount<<endl;
                //contour -> set(x,y,0);
                traceContour(x, y, contour);
            }
    
    OutFile2<<"Delta encoding contour pixels needs: "<<ByteCount<<" Bytes."<<endl;
   
    OutFile2.close();
}

void Image::RemoveRepeatPoints()
{
    ifstream ifs("controlPoint.txt"); 
    string str;
    ofstream CPandIndex, Indexstream, IndexSize;
    CPandIndex.open("CPandIndex.txt");
    Indexstream.open("Index.txt");
    IndexSize.open("../output.txt", ios_base::app);
    int CPnum, degree, SampleNum;
    vector<Vector3<int>> storeCPs;
    Vector3<int> currentPoint;
    vector<int> Index;
    int index = 0;
    bool find = false;

    ifs >> str;
    int BranchSize = (int)(atof(str.c_str()));
    BranchSize += SplitNum;
    CPandIndex << BranchSize <<endl;    
    while(BranchSize--)
    {
        ifs >> str;
        CPnum = (int)(atof(str.c_str()));
        
        ifs >> str;
        degree = (int)(atof(str.c_str()));
        
        ifs >> str;
        SampleNum = (int)(atof(str.c_str()));
        CPandIndex << CPnum <<" "<< degree <<" "<< SampleNum <<" ";

        while(CPnum--)  
        {
            ifs >> str;
            currentPoint[0] = (int)(atof(str.c_str()));
            
            ifs >> str;
            currentPoint[1] = (int)(atof(str.c_str()));
        
            ifs >> str;
            currentPoint[2] = (int)(atof(str.c_str()));

             
           //find//
           if(storeCPs.empty())//first time
           {
                CPandIndex<<currentPoint[0]<<" "<<currentPoint[1]<<" "<<currentPoint[2]<<" ";
                index++;
                storeCPs.push_back(currentPoint);
           }
           else
           {
               for(auto it = storeCPs.begin();it!=storeCPs.end();it++)
               {
                   Vector3<int> point = *it;
                   if(abs(currentPoint[0]-point[0]) < 2 && abs(currentPoint[1]-point[1]) < 2 && abs(currentPoint[2]-point[2]) < 2)
                    {
                        Index.push_back(BranchSize);
                        Index.push_back(index);
                        int id = distance(storeCPs.begin(), it);
                        Index.push_back(id);
                        find = true;
                        break;
                    }
               }
                if(!find)
                {
                    CPandIndex<<currentPoint[0]<<" "<<currentPoint[1]<<" "<<currentPoint[2]<<" ";
                    index++;
                    storeCPs.push_back(currentPoint);
                }
                find = false;  
           }
           ////end
        }
        index = 0;
        CPandIndex<<endl;
    }
    ifs.close();
    for(auto it = Index.begin();it!=Index.end();it++)
        Indexstream << *it <<" ";
    
    IndexSize<<"IndexSize: "<<Index.size()/3<<endl;
    IndexSize.close();
    CPandIndex.close();
    Indexstream.close();
    
}


/**
* Calculate the skeleton branches for binary image.
*/
skel_tree_t* Image::computeBinarySkel()
{
    int fboSize = initialize_skeletonization(im);
    FIELD<float>* imDt = im->dupe();
    clock_t start = clock();
    skelImp = computeSkeleton(imDt);
    clock_t endtime = clock();
    cout<<"run time: "<<(double)(endtime-start)/CLOCKS_PER_SEC<<" s" <<endl;

    int max = 0;
    //NewRemoveDoubleSkel();//might be broke in junction point.
    
    FIELD<float>* skel = new FIELD<float>(skelImp->dimX(), skelImp->dimY());
    for (int i = 0; i < skelImp->dimX(); ++i) {
        for (int j = 0; j < skelImp->dimY(); ++j) {
            bool is_skel_point = skelImp->value(i,j);
            skel->set(i, j, is_skel_point ? 255 : 0);
            if (skelImp->value(i,j) > max) max = skelImp->value(i,j);
        }
    } 

    removeDoubleSkel(skel);
    removeRtZero(skel,imDt);
    
    analyze_cca(skel);//delete short branches further.
    
    skel -> writePGM("origskel.pgm");
    //count
    int count = 0;
    for (int y = 0; y < skel->dimY(); ++y) 
        for (int x = 0; x < skel->dimX(); ++x) 
            if (skel->value(x, y))
				count++;
    cout<< "MAT count: "<<count<<endl;
    //ofstream OutFile1;
    OutFile2.open("../output.txt",ios_base::app);
    OutFile2 << "Storing MAT pixels needs: " << 6 * count + 4 <<" Bytes."<<endl;         
    
    //imDt -> writePGM("dt.pgm");
    FIELD<float>* newSkel = new FIELD<float> (*skel);
    
    //segment = new FIELD<float>(skel->dimX(), skel->dimY());
    //removeDoubleSkel(newSkel);
    int seq = 0;
   /////segment and store into the SampleSets;
    diagonal  = sqrt((float)(newSkel->dimX() * newSkel->dimX() + newSkel->dimY() * newSkel->dimY()));
    
	for (int y = 0; y < newSkel->dimY(); ++y) {
        for (int x = 0; x < newSkel->dimX(); ++x) {
            if (newSkel->value(x, y) > 0) {
                vector<Vector3<float>> Branch;
                tracePath(x, y, newSkel, imDt, Branch, seq);//passed by reference
                seq++;
            }
        }
    }
    
    OutFile2.close();    
    
    ////fit with spline///
    cout<<"BranchSet.size()--"<<BranchSet.size()<<endl;
    ///BSplineCurveFitterWindow3 spline(BranchSet, hausdorff, diagonal);
    BSplineCurveFitterWindow3 spline;
    SplitNum = spline.SplineFit(BranchSet, hausdorff, diagonal, connection);
    // Trace the tree and store it.
    skel_tree_t* tree = traceLayer(skel, imDt); 

    return tree;
}
