/*
 * File:   ImageDecoder.cpp
 * Author: yuri
 *
 * Created on June 13, 2011, 1:32 PM
 */
 
#include <vector>
#include "./include/io.hpp"
#include "./include/messages.h"
#include "./include/DataManager.hpp"
#include "./include/ImageDecoder.hpp"
#include "../shared/FastAC/arithmetic_codec.h"
#include "fileio/fileio.hpp"
#include <assert.h>
#include <fstream>
#include <bitset>
#include <climits>
#include <iostream>
#include <cmath>
#include "SplineGenerate/BSplineCurveFitter/BSplineCurveFitterWindow3.h"
//#include <string>
using namespace std;

ImageDecoder::ImageDecoder()
{
    width = -1;
    height = -1;
}

ImageDecoder::ImageDecoder(const ImageDecoder& orig) {
}

int MAXR = 0;

ImageDecoder::~ImageDecoder() {
}

/* Add a point to the current datatype. This is used for both starting points, and neighbouring points.
 * Identify a starting point by setting "point" to 0. This cannot happen for neighbouring points, as that
 * implicates the radius is 0.*/
void ImageDecoder::addPoint(char point, int &x, int &y, int &r, path_t *path) {
    /* If we are a neighbouring point: decode the character and walk in the direction. Modify the original
     * x,y,r values. */
    // Layout of point is 0xxyyrrr where xx are the bits denoting delta_x [0, 1, 2]
    // yy denotes delta_y [0, 1, 2] and rrr denote delta_radius [0, 1, 2, 3, 4]
    int8_t dx = ((point >> 5) & 0x3) - 1;
    int8_t dy = ((point >> 3) & 0x3) - 1;
    int8_t dr = (point & 0x7) - 2;
    x += dx;
    y += dy;
    r += dr;


    // if (r < 0) { cerr << "Invalid radius size, must be >0 [r=" << r << "]" << endl; exit(-1);}

    /* Add to datatype */
    path->push_back(coord3D_t(x, y, r));

}

const char* ImageDecoder::get_comp_method_string(COMPRESS_MODE mode) {
    switch (mode) {
    case COMPRESS_ZLIB:
        return "zlib";
        break;
    case COMPRESS_LZHAM:
        return "lzham";
        break;
    case COMPRESS_BSC:
        return "bsc";
        break;
    case COMPRESS_CSC:
        return "csc";
        break;
    case COMPRESS_LZMA:
        return "xz";
        break;
    case COMPRESS_LZMA2:
        return "lzma2";
        break;
    case COMPRESS_BROTLI:
        return "brotli";
        break;
    case COMPRESS_ZPAQ:
        return "zpaq";
        break;
    case COMPRESS_BZIP2:
        return "bzip2";
        break;
    case COMPRESS_ZSTD:
        return "zstd";
        break;
    default:
        // Uses some unsupported method or file is mangled.
        return NULL;
        break;
    }
}

layer_t* ImageDecoder::decode_layer()
{
    layer_t *tree = new layer_t();
    uint16_t numPaths = READUINT16(dat_it);
    uint8_t current;
    if (numPaths == 0xFFFF)
        cout<<"fail numPath!!"<<endl;

    of_uncompressed << "Num children " << +numPaths << endl;
    
    for (int i = 0; i < numPaths; ++i) {
     
        path_t* path = new path_t();
        /* Read first -full- point */
        /* Seriously, y tho*/
        uint8_t x1 = READUINT8(dat_it);
        uint8_t x2 = READUINT8(dat_it);
        int x = (x1 << 8) + x2;
        uint8_t y1 = READUINT8(dat_it);
        uint8_t y2 = READUINT8(dat_it);
        int y = (y1 << 8) + y2;
        uint8_t r1 = READUINT8(dat_it);
        uint8_t r2 = READUINT8(dat_it);
        int r = (r1 << 8) + r2;
        path->push_back(coord3D_t(x, y, r));
      
        bool end = false;
        while (!end) {
            of_uncompressed << x << " - " << y << " - " << r << endl;
            current = READUINT8(dat_it);
            if (current == END_TAG) {
                of_uncompressed << "End" << endl;
                end = true;
            }

            /* End of branch? */
            else if (current == FORK_TAG) {
                uint8_t goBack_1 = READUINT8(dat_it);
                uint8_t goBack_2 = READUINT8(dat_it);
                uint16_t goBack = (goBack_1 << 8) + goBack_2; // TODO(maarten): why is this necessary?
                of_uncompressed << "Fork - " << (goBack) << endl;
                for (unsigned int q = 0; q < goBack; ++q) {
                    tree->push_back((path->back()));
                    path->pop_back();
                }
                x = path->back().first;
                y = path->back().second;
                r = path->back().third;
                
            } else {
                if ((current < 128)) {
                    addPoint(current, x, y, r, path);
                } else {
                    if (current == 128) {
                        if (bits_dx == 0) {
                            int8_t dx_8 = READINT8(dat_it);
                            x += dx_8;
                        } else {
                            uint8_t dx_1_16 = READUINT8(dat_it);
                            uint8_t dx_2_16 = READUINT8(dat_it);
                            x += ((dx_1_16 << 8) + dx_2_16);
                            if (x > width)
                                x -= (1 << 16);
                        }
                        if (bits_dy == 0) {
                            int8_t dy_8 = READINT8(dat_it);
                            y += dy_8;
                        } else {
                            uint8_t dy_1_16 = READUINT8(dat_it);
                            uint8_t dy_2_16 = READUINT8(dat_it);
                            y += ((dy_1_16 << 8) + dy_2_16);
                            if (y > width)
                                y -= (1 << 16);

                        }
                        if (bits_dr == 0) {
                            int8_t dr_8 = READINT8(dat_it);
                            r += dr_8;
                        } else {
                            uint8_t dr_1_16 = READUINT8(dat_it);
                            uint8_t dr_2_16 = READUINT8(dat_it);
                            r += ((dr_1_16 << 8) + dr_2_16);
                            if (r > width)
                                r -= (1 << 16);

                        }
                    } else {
                        uint16_t num = (current << 8) + READUINT8(dat_it);
                        int8_t dr = ((num & 0x1F)) - 15;
                        int8_t dy = (((num >> 5) & 0x1F)) - 15;
                        int8_t dx = (((num >> 10) & 0x1F)) - 15;
                        x += dx;
                        y += dy;
                        r += dr;
                    }

                    if (r < 0) { cerr << "Invalid radius size, must be >0 [r=" << r << "]" << endl; exit(-1);}

                    /* Add to datatype */
                    path->push_back(coord3D_t(x, y, r));
                }
            }

        }

        /* Copy path to the layer of the image */
        for (unsigned int i = 0; i < path->size(); ++i) {
           /* if (intensity ==249)
            {
                if(i == 0){     //just for debug
                PRINT(MSG_NORMAL, "pathx %d\n", (*path)[1].first);
                PRINT(MSG_NORMAL, "pathy %d\n", (*path)[1].second);
                PRINT(MSG_NORMAL, "pathz %d\n", (*path)[1].third);
                }
            } */
                
            tree->push_back((*path)[i]);
        }

    }
    return tree;
}

bool ImageDecoder::decode_layer(image_t** image_ref, int intensity) {
    // uint8_t x1 = READUINT8(dat_it);
    // uint8_t x2 = READUINT8(dat_it);
   if (intensity == 255) intensity = 254;//Wang. I don't know why when the intensity equals to 255, it doesn't work, so here I just change it to 254 since it isn't distinguishable.
    image_t* image = *image_ref;
    path_t* layer = (*image)[intensity];
    // uint16_t numPaths =  (x1 << 8) + x2;
    uint16_t numPaths = READUINT16(dat_it);
    uint8_t current;
    if (numPaths == 0xFFFF)
        return false;

    of_uncompressed << "Intensity " << +intensity << " - " << "Num children " << +numPaths << endl;
    
    for (int i = 0; i < numPaths; ++i) {
     
        path_t* path = new path_t();
        /* Read first -full- point */
        /* Seriously, y tho*/
        uint8_t x1 = READUINT8(dat_it);
        uint8_t x2 = READUINT8(dat_it);
        int x = (x1 << 8) + x2;
        uint8_t y1 = READUINT8(dat_it);
        uint8_t y2 = READUINT8(dat_it);
        int y = (y1 << 8) + y2;
        uint8_t r1 = READUINT8(dat_it);
        uint8_t r2 = READUINT8(dat_it);
        int r = (r1 << 8) + r2;
        path->push_back(coord3D_t(x, y, r));
      
        bool end = false;
        while (!end) {
            of_uncompressed << x << " - " << y << " - " << r << endl;
            current = READUINT8(dat_it);
            if (current == END_TAG) {
                of_uncompressed << "End" << endl;
                end = true;
            }

            /* End of branch? */
            else if (current == FORK_TAG) {
                uint8_t goBack_1 = READUINT8(dat_it);
                uint8_t goBack_2 = READUINT8(dat_it);
                uint16_t goBack = (goBack_1 << 8) + goBack_2; // TODO(maarten): why is this necessary?
                of_uncompressed << "Fork - " << (goBack) << endl;
                for (unsigned int q = 0; q < goBack; ++q) {
                    layer->push_back((path->back()));
                    path->pop_back();
                }
                x = path->back().first;
                y = path->back().second;
                r = path->back().third;
            } else {
                if ((current < 128)) {
                    addPoint(current, x, y, r, path);
                } else {
                    if (current == 128) {
                        if (bits_dx == 0) {
                            int8_t dx_8 = READINT8(dat_it);
                            x += dx_8;
                        } else {
                            uint8_t dx_1_16 = READUINT8(dat_it);
                            uint8_t dx_2_16 = READUINT8(dat_it);
                            x += ((dx_1_16 << 8) + dx_2_16);
                            if (x > width)
                                x -= (1 << 16);
                        }
                        if (bits_dy == 0) {
                            int8_t dy_8 = READINT8(dat_it);
                            y += dy_8;
                        } else {
                            uint8_t dy_1_16 = READUINT8(dat_it);
                            uint8_t dy_2_16 = READUINT8(dat_it);
                            y += ((dy_1_16 << 8) + dy_2_16);
                            if (y > width)
                                y -= (1 << 16);

                        }
                        if (bits_dr == 0) {
                            int8_t dr_8 = READINT8(dat_it);
                            r += dr_8;
                        } else {
                            uint8_t dr_1_16 = READUINT8(dat_it);
                            uint8_t dr_2_16 = READUINT8(dat_it);
                            r += ((dr_1_16 << 8) + dr_2_16);
                            if (r > width)
                                r -= (1 << 16);

                        }
                    } else {
                        uint16_t num = (current << 8) + READUINT8(dat_it);
                        int8_t dr = ((num & 0x1F)) - 15;
                        int8_t dy = (((num >> 5) & 0x1F)) - 15;
                        int8_t dx = (((num >> 10) & 0x1F)) - 15;
                        x += dx;
                        y += dy;
                        r += dr;
                    }

                    if (r < 0) { cerr << "Invalid radius size, must be >0 [r=" << r << "]" << endl; exit(-1);}

                    /* Add to datatype */
                    path->push_back(coord3D_t(x, y, r));
                }
            }

        }
        
        /* Copy path to the layer of the image */
        for (unsigned int i = 0; i < path->size(); ++i) {
           /* if (intensity ==249)
            {
                if(i == 0){     //just for debug
                PRINT(MSG_NORMAL, "pathx %d\n", (*path)[1].first);
                PRINT(MSG_NORMAL, "pathy %d\n", (*path)[1].second);
                PRINT(MSG_NORMAL, "pathz %d\n", (*path)[1].third);
                }
            } */
                
            layer->push_back((*path)[i]);
        }

    }
    //PRINT(MSG_NORMAL, "33\n");
    return true;
}

image_t* ImageDecoder::decode_plane(vector<int>& levels) {
    image_t* img = new image_t();
   
    levels.clear();
    int num_levels = READUINT32(dat_it);
    for (int i = 0; i < num_levels; ++i) {
        int level = READUINT8(dat_it);
        if (i == 0)
            levels.push_back(level);
        else
            levels.push_back(levels.at(i - 1) + level);
    }
    for (int i = 0; i < 0xff; ++i) {
        path_t *p = new path_t();
        img->push_back(p);
    }
    of_uncompressed.open("uncompressed.sir", ios_base::out | ios::binary);
    //PRINT(MSG_NORMAL, "1\n");
    for (auto intensity : levels) {
        decode_layer(&img, intensity);
    }
    of_uncompressed.close();
    return img;
}

void ImageDecoder::loadSMAT_R(const char *fname){

    unsigned int nEl = readFile<unsigned char>(fname, &data);
    if (nEl == 0) {
        PRINT(MSG_ERROR, "Could not open file for Image Decoder.\n");
        exit(1);
    }

    int UseSquash = (int)READUINT8(data);//First Byte store the sign that if we use squash compression.
    
    // Initialize iterator 
    dat_it = data;

    ofstream OutFile;
    OutFile.open("ControlPoints.txt");//wang
    int x,y,dt;
    int BranchSize;
    uint8_t x1,x2;
    vector<int> Index;
    vector<Vector3<int>> storeCPs;
    Vector3<int> currentPoint, repeatPoint;


    x1 = READUINT8(dat_it);
    x2 = READUINT8(dat_it);
    x = ((x1 << 8) + x2);//width
    x1 = READUINT8(dat_it);
    x2 = READUINT8(dat_it);
    y = ((x1 << 8) + x2);//height
    
    OutFile<<x<<" "<<y<<endl;

    x1 = READUINT8(dat_it);
    BranchSize = (int) x1 ;//BranchSize

    while (true)
    {
        x1 = READUINT8(dat_it);
        x = (int) x1;
        
        if(x == 255) break;
        else
        {
            Index.push_back(x);
            x1 = READUINT8(dat_it);
            x2 = READUINT8(dat_it);
            int second = (int) ((x1 >> 4) & 0xF);//first 4 bits
            Index.push_back(second);
            int third = ((x1 & 0xF) << 8) + x2;
            Index.push_back(third);
        }  
    }

    while (BranchSize--)
    {
        uint8_t firstB = READUINT8(dat_it);
        int degree = (int) (firstB  & 0xF);
        int NumCP = (int) ((firstB >> 4) & 0xF);
        
        OutFile<<NumCP<<" "<<degree<<" ";

        x1 = READUINT8(dat_it);
        x2 = READUINT8(dat_it);
        x = (x1 << 8) + x2;
        OutFile<<x<<" ";//numSample

        ////judge if this branch has repeat points

    while(true){

        if(Index.empty())  break;
        auto it = Index.begin();
        if(BranchSize == *it)//this branch has repeat points.
        {
            NumCP--;
            Index.erase(Index.begin());
            it = Index.begin();
            int before = *it;
            while(before--)
            {
                NumCP--;
                x1 = READUINT8(dat_it);
                x2 = READUINT8(dat_it);
                x = ((x1 << 8) + x2);
                currentPoint[0] = x > 32768 ? (x-65536) : x;
                x1 = READUINT8(dat_it);
                x2 = READUINT8(dat_it);
                y = (x1 << 8) + x2;
                currentPoint[1] = y > 32768 ? (y-65536) : y;
                x1 = READUINT8(dat_it);
                x2 = READUINT8(dat_it);
                dt = (x1 << 8) + x2;
                currentPoint[2] = dt > 32768 ? (dt-65536) : dt;
                OutFile<<currentPoint[0]<<" "<<currentPoint[1]<<" "<<currentPoint[2]<<" ";
                storeCPs.push_back(currentPoint);
            }
            Index.erase(Index.begin());
            it = Index.begin();
            repeatPoint = storeCPs.at(*it);
            OutFile<<repeatPoint[0]<<" "<<repeatPoint[1]<<" "<<repeatPoint[2]<<" ";
            Index.erase(Index.begin());
            
        }
        else break;
    }

        while(NumCP--)
        {
            x1 = READUINT8(dat_it);
            x2 = READUINT8(dat_it);
            x = ((x1 << 8) + x2);
            currentPoint[0] = x > 32768 ? (x-65536) : x;
            x1 = READUINT8(dat_it);
            x2 = READUINT8(dat_it);
            y = (x1 << 8) + x2;
            currentPoint[1] = y > 32768 ? (y-65536) : y;
            x1 = READUINT8(dat_it);
            x2 = READUINT8(dat_it);
            dt = (x1 << 8) + x2;
            currentPoint[2] = dt > 32768 ? (dt-65536) : dt;
            OutFile<<currentPoint[0]<<" "<<currentPoint[1]<<" "<<currentPoint[2]<<" ";
            storeCPs.push_back(currentPoint);
        }
        OutFile<<endl;
    }
}


void ImageDecoder::loadSMAT_(const char *fname){

    unsigned int nEl = readFile<unsigned char>(fname, &data);
    if (nEl == 0) {
        PRINT(MSG_ERROR, "Could not open file for Image Decoder.\n");
        exit(1);
    }

    int UseSquash = (int)READUINT8(data);//First Byte store the sign that if we use squash compression.
    
    // Initialize iterator 
    dat_it = data;

    ofstream OutFile;
    OutFile.open("ControlPoints.txt");//wang
    int x,y,dt;
    int last_x = 0, last_y = 0, last_dt = 0;
    uint8_t x1,x2;

    x1 = READUINT8(dat_it);
    x2 = READUINT8(dat_it);
    x = ((x1 << 8) + x2);//width
    x1 = READUINT8(dat_it);
    x2 = READUINT8(dat_it);
    y = ((x1 << 8) + x2);//height
    
    OutFile<<x<<" "<<y<<endl;
    while (true)
    {
        uint8_t firstB = READUINT8(dat_it);
        //cout<<firstB<<endl;
        if ((int)firstB == 255) break;
        int degree = (int) (firstB  & 0xF);
        int NumCP = (int) ((firstB >> 4) & 0xF);
        
        OutFile<<NumCP<<" "<<degree<<" ";

        x1 = READUINT8(dat_it);
        x2 = READUINT8(dat_it);
        x = (x1 << 8) + x2;
        OutFile<<x<<" ";//numSample
        while(NumCP--)
        {
            x1 = READUINT8(dat_it);
            x2 = READUINT8(dat_it);
            x = ((x1 << 8) + x2);
            //cout<<"x----------"<<x<<endl;
            if (x == 25700)
            {
                x1 = READUINT8(dat_it);
                x2 = READUINT8(dat_it);
                x = ((x1 << 8) + x2);
                x = x > 32768 ? (x-65536) : x;
                x += last_x;

                x1 = READUINT8(dat_it);
                x2 = READUINT8(dat_it);
                y = (x1 << 8) + x2;
                y = y > 32768 ? (y-65536) : y;
                y += last_y;

                x1 = READUINT8(dat_it);
                x2 = READUINT8(dat_it);
                dt = (x1 << 8) + x2;
                dt = dt > 32768 ? (dt-65536) : dt;
                dt += last_dt;
            }
            else
            {
                x = (int)x1;
                x = x > 128 ? (x-256) : x;
                x += last_x;

                y = (int)x2;
                y = y > 128 ? (y-256) : y;
                y += last_y;

                x1 = READUINT8(dat_it);
                dt = (int)x1;
                dt = dt > 128 ? (dt-256) : dt;
                dt += last_dt;
            }
            last_x = x;
            last_y = y;
            last_dt = dt;
            OutFile<<x<<" "<<y<<" "<<dt<<" ";
        }
        OutFile<<endl;
    }
}

void ImageDecoder::loadSMAT(const char *fname){

    unsigned int nEl = readFile<unsigned char>(fname, &data);
    if (nEl == 0) {
        PRINT(MSG_ERROR, "Could not open file for Image Decoder.\n");
        exit(1);
    }

    int UseSquash = (int)READUINT8(data);//First Byte store the sign that if we use squash compression.
    
    // Initialize iterator 
    dat_it = data;

    ofstream OutFile;
    OutFile.open("ControlPoints.txt");//wang
    int x,y,dt;
    uint8_t x1,x2;

    x1 = READUINT8(dat_it);
    x2 = READUINT8(dat_it);
    x = ((x1 << 8) + x2);//width
    x1 = READUINT8(dat_it);
    x2 = READUINT8(dat_it);
    y = ((x1 << 8) + x2);//height
    
    OutFile<<x<<" "<<y<<endl;
    while (true)
    {
        uint8_t firstB = READUINT8(dat_it);
        //cout<<firstB<<endl;
        if ((int)firstB == 255) break;
        int degree = (int) (firstB  & 0xF);
        int NumCP = (int) ((firstB >> 4) & 0xF);
        
        OutFile<<NumCP<<" "<<degree<<" ";

        x1 = READUINT8(dat_it);
        x2 = READUINT8(dat_it);
        x = (x1 << 8) + x2;
        OutFile<<x<<" ";//numSample
        while(NumCP--)
        {
            x1 = READUINT8(dat_it);
            x2 = READUINT8(dat_it);
            x = ((x1 << 8) + x2);
            x = x > 32768 ? (x-65536) : x;
            x1 = READUINT8(dat_it);
            x2 = READUINT8(dat_it);
            y = (x1 << 8) + x2;
            y = y > 32768 ? (y-65536) : y;
            x1 = READUINT8(dat_it);
            x2 = READUINT8(dat_it);
            dt = (x1 << 8) + x2;
            dt = dt > 32768 ? (dt-65536) : dt;
            OutFile<<x<<" "<<y<<" "<<dt<<" ";
        }
        OutFile<<endl;
    }
}

FIELD<float>* ImageDecoder::loadSample(int SuperResolution) {////
    
    //COLORSPACE colorspace = COLORSPACE::GRAY;
    clear_color = 0;
    binary = new layer_t();
    ifstream ifs("sample.txt"); 
    string str;
    ifs >> str;
    width = (int)atof(str.c_str());
    cout<<"width: "<<width<<endl;
    ifs >> str;
    height = (int)atof(str.c_str());        
    int x,y,dt,imp;
    FIELD<float>* impmap = new FIELD<float>(width,height);
    
    while(ifs)
    { 
        ifs >> str;
        x = (round)(atof(str.c_str())*SuperResolution);//sub-pixel
        ifs >> str;
        y = (round)(atof(str.c_str())*SuperResolution);
        ifs >> str;
        dt = (round)(atof(str.c_str())*SuperResolution);
        //ifs >> str;//4D
        //imp = (round)(atof(str.c_str()));//4D

        if(ifs.fail())  break;
        //impmap->set(x,y,imp);//4D
        binary->push_back(coord3D_t(x, y, dt));

    }   
        //impmap->writePGM("impmap.pgm");//4D
        ifs.close();
    return impmap;
}
