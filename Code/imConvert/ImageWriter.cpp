/*
 * File:   ImageWriter.cpp
 * Author: yuri
 *
 * Created on June 6, 2011, 2:59 PM
 */
#include <omp.h>
#include "include/ImageWriter.hpp"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <bitset>
#include <math.h>
#include "include/messages.h"
#include "include/io.hpp"
#include "include/SkelNode.hpp"
#include "fileio/fileio.hpp"
#include "include/HuffNode.hpp"
#include "bcl/rle.h"
#include <algorithm>
#include <set>
#include <queue>
#include <boost/dynamic_bitset.hpp>
#include "include/ImageEncoder.hpp"
#include "SplineGenerate/BSplineCurveFitter/BSplineCurveFitterWindow3.h"

using namespace std;

enum IW_STATE { INVALID = -1,
                FILE_OPEN,
                HEADER_WRITTEN,
                FILE_SAVED,
                FILE_CLOSE
              } state;

string COMP_METHOD_NAME;
string ENCODING;
int COMP_LEVEL;
//ofstream fout;

bool icompare_pred(unsigned char a, unsigned char b) {
    return std::tolower(a) == std::tolower(b);
}

bool iequals(std::string const& a, std::string const& b) {
    if (a.length() == b.length()) {
        return std::equal(b.begin(), b.end(),
                          a.begin(), icompare_pred);
    } else {
        return false;
    }
}
ImageWriter::ImageWriter(const char *fname) {
    of.open(fname, ios_base::out | ios::binary);
    of_uncompressed.open("uncompressed.sir", ios_base::out | ios::binary);

    if (of.good() && of_uncompressed.good())
    { state = FILE_OPEN; }
    else {
        PRINT(MSG_ERROR, "Could not open file for output\n");
        exit(1);
    }
    if (iequals(COMP_METHOD_NAME, "lzma")) {
        comp_mode = COMPRESS_LZMA;
    } else if (iequals(COMP_METHOD_NAME, "csc")) {
        comp_mode = COMPRESS_CSC;
    } else if (iequals(COMP_METHOD_NAME, "lzham")) {
        comp_mode = COMPRESS_LZHAM;
    } else if (iequals(COMP_METHOD_NAME, "brotli")) {
        comp_mode = COMPRESS_BROTLI;
    } else if (iequals(COMP_METHOD_NAME, "zpaq")) {
        comp_mode = COMPRESS_ZPAQ;
    } else if (iequals(COMP_METHOD_NAME, "bzip2")) {
        comp_mode = COMPRESS_BZIP2;
    } else if (iequals(COMP_METHOD_NAME, "lzma2")) {
        comp_mode = COMPRESS_LZMA2;
    } else if (iequals(COMP_METHOD_NAME, "bsc")) {
        comp_mode = COMPRESS_BSC;
    } else if (iequals(COMP_METHOD_NAME, "zlib")) {
        comp_mode = COMPRESS_ZLIB;
    } else if (iequals(COMP_METHOD_NAME, "zstd")) {
        comp_mode = COMPRESS_ZSTD;
    } else {
        PRINT(MSG_NORMAL, "Compression method not recognized. Defaulting to LZMA\n");
        comp_mode = COMPRESS_LZMA;
    }
    if (iequals(ENCODING, "huffman")) {
        encode_mode = (ENCODE_MODE::HUFFMAN);
    } else if (iequals(ENCODING, "canonical")) {
        encode_mode = (ENCODE_MODE::CANONICAL);
    } else if (iequals(ENCODING, "unitary")) {
        encode_mode = (ENCODE_MODE::UNITARY);
    } else if (iequals(ENCODING, "exp_goulomb")) {
        encode_mode = (ENCODE_MODE::EXP_GOULOMB);
    } else if (iequals(ENCODING, "arithmetic")) {
        encode_mode = (ENCODE_MODE::ARITHMETIC);
    } else if (iequals(ENCODING, "predictive")) {
        encode_mode = (ENCODE_MODE::PREDICTIVE);
    } else if (iequals(ENCODING, "compact")) {
        encode_mode = (ENCODE_MODE::COMPACT);
    } else if (iequals(ENCODING, "raw")) {
        encode_mode = (ENCODE_MODE::RAW);
    } else if (iequals(ENCODING, "mtf")) {
        encode_mode = (ENCODE_MODE::MTF);
    } else {
        encode_mode = (ENCODE_MODE::TRADITIONAL);
    }
}

string ImageWriter::get_compression_name(COMPRESS_MODE mode) {
    string res = "";
    switch (mode) {
    case COMPRESS_LZMA:
        res = "xz";
        break;
    case COMPRESS_LZMA2:
        res = "lzma2";
        break;
    case COMPRESS_BROTLI:
        res = "brotli";
        break;
    case COMPRESS_ZPAQ:
        res = "zpaq";
        break;
    case COMPRESS_BSC:
        res = "bsc";
        break;
    case COMPRESS_BZIP2:
        res = "bzip2";
        break;
    case COMPRESS_ZLIB:
        res = "zlib";
        break;
    case COMPRESS_LZHAM:
        res = "lzham";
        break;
    case COMPRESS_CSC:
        res = "csc";
        break;
    case COMPRESS_ZSTD:
        res = "zstd";
        break;
    default:
        res = "lzma";
        break;
    }
    return res;
}

void ImageWriter::save() {
    ret_struct val_struct;
    val_struct.mode = comp_mode;
    val_struct.compress_algo = get_compression_name(comp_mode);
    PRINT(MSG_VERBOSE, "Compressing using %s\n", val_struct.compress_algo.c_str());
    val_struct.in_data = ofBuffer.str();
    val_struct.in_size = ofBuffer.str().size();
   
   //Always do not use squash library.
    of.write((const char *) "00000000", sizeof(uint8_t));
    of_uncompressed << "squash?" << "00000000" << '\n';
    auto out = reinterpret_cast<unsigned char*>(const_cast<char*>(ofBuffer.str().c_str()));
    PRINT(MSG_NORMAL, "NoSquash: Bits / pixel: %3.4f\n", ofBuffer.str().size() / (double)(num_pixels));
    of.write((const char *) &out[0], ofBuffer.str().size());
    
    state = FILE_SAVED;

    of.close();
    of_uncompressed.close();
}

ImageWriter::~ImageWriter() {
    if (state != FILE_SAVED) {
        save();
    }
    of.close();
    of_uncompressed.close();
}

void ImageWriter::writeBits(ostream& os, const string& str) {
    for (unsigned i = 0; i < str.length(); i += 8) {
        string sub = str.substr(i, 8);
        uint8_t bit_rep = bitset<8>(sub).to_ulong();
        os.write((const char*)(&bit_rep), sizeof(char));
    }
}
 
   /* string encode_start_point(int x, int y, int r, int imp) {
        return bitset<16>(x).to_string() + bitset<16>(y).to_string() + bitset<16>(r).to_string() + bitset<16>(imp).to_string();
    }////4D*/
   string encode_start_point(int x, int y, int r) {
        return bitset<16>(x).to_string() + bitset<16>(y).to_string() + bitset<16>(r).to_string();
    }////3D
    string encode_8bit_point(int x, int y, int r) {
        return bitset<8>(x).to_string() + bitset<8>(y).to_string() + bitset<8>(r).to_string();
    }////3D
    string encode_16bit_point(int x, int y, int r) {
        return bitset<16>(x).to_string() + bitset<16>(y).to_string() + bitset<16>(r).to_string();
    }////3D

void ImageWriter::writeCPandID(unsigned int width, unsigned int height) {
    
    ifstream ifs("CPandIndex.txt"); 
    string str;
    int x,y,dt,imp; ////if we use uint16_t, then we cannot represent negative value.(-5 turns to 65531)
    int CPnum, degree;
    
    num_pixels = width * height;

    of_uncompressed << width << " - " << height << '\n';
    writeBits(ofBuffer, bitset<16>(width).to_string() + bitset<16>(height).to_string());

    ifs >> str;
    int BranchSize = (int)(atof(str.c_str()));
    of_uncompressed << BranchSize << '\n';
    writeBits(ofBuffer, bitset<8>(BranchSize).to_string());

    if( BranchSize > 255)
    {
        cout<<"ERROR!!!!Branchsize is bigger than 256, Need to redistribute!"<<endl;
        return;
    }
///write index
    ifstream ifs1("Index.txt");
    string str1;
    int index;
    vector<int> Index;
    string firstFourBits,lastFourBits;

    while(ifs1)
    {
        ifs1 >> str1;
        if(ifs1.fail()) break;

        index = (int)(atof(str1.c_str()));///
        of_uncompressed << index << '\n';
        writeBits(ofBuffer, bitset<8>(index).to_string());
        Index.push_back(index);
        
        //---------
        ifs1 >> str1;
        index = (int)(atof(str1.c_str()));///
        Index.push_back(index);
        firstFourBits = bitset<4>(index).to_string();
        of_uncompressed << index << " - ";
        //---------
        ifs1 >> str1;
        index = (int)(atof(str1.c_str()));///
        of_uncompressed << index << '\n';
        Index.push_back(index);
        lastFourBits = bitset<12>(index).to_string();
        writeBits(ofBuffer, firstFourBits + lastFourBits);

    }
     of_uncompressed << "11111111" << '\n';
     writeBits(ofBuffer, "11111111");
///write branches
    while(BranchSize--)
    {
        ifs >> str;
        CPnum = (int)(atof(str.c_str()));
        firstFourBits = bitset<4>(CPnum).to_string();
    
        ifs >> str;
        degree = (int)(atof(str.c_str()));
        lastFourBits = bitset<4>(degree).to_string();

        //of_uncompressed << firstFourBits << " - " << lastFourBits << '\n';
        of_uncompressed << CPnum << " - " << degree << '\n';
        writeBits(ofBuffer, firstFourBits + lastFourBits);

        ifs >> str;
        string SampleNum = bitset<16>((int)(atof(str.c_str()))).to_string();
        of_uncompressed << (int)(atof(str.c_str())) << '\n';
        writeBits(ofBuffer, SampleNum);

////judge if this branch has repeat points
        if(!Index.empty()) {
            while(true){
                auto it = Index.begin();
                if(BranchSize == *it)//this branch has repeat points.
                {
                    CPnum--;
                    Index.erase(Index.begin());
                    Index.erase(Index.begin());
                    Index.erase(Index.begin());
                    if(Index.empty())  break;
                }
                else break;
            }
    } 
        while(CPnum--)  
        {
            
            ifs >> str;
            x = (atof(str.c_str()));
            ifs >> str;
            y = (atof(str.c_str()));
            ifs >> str;
            dt = (atof(str.c_str()));
            //ifs >> str;////4D
            //imp = (atof(str.c_str()));

            of_uncompressed << x << " - " << y << " - " << dt << '\n';
            writeBits(ofBuffer, encode_16bit_point(x,y,dt));
        } 
        
    }  
}
//improve method- 8-bit version.
void ImageWriter::writeCP_(unsigned int width, unsigned int height) {
    ofstream  AddToOutput;
    AddToOutput.open("../output.txt", ios_base::app);

    ifstream ifs("controlPoint.txt"); 
    string str;
    int imp; ////if we use uint16_t, then we cannot represent negative value.(-5 turns to 65531)
    int CPnum, degree;
    int last_x = 0, last_y = 0, last_dt = 0;
    int encode_x, encode_y, encode_dt;
    
    vector<Vector3<int>> storeCPs;
    Vector3<int> currentPoint;
    int RepeatNum = 0; bool find = false;

    num_pixels = width * height;

    of_uncompressed << width << " - " << height << '\n';
    writeBits(ofBuffer, bitset<16>(width).to_string() + bitset<16>(height).to_string());

    ifs >> str;//branchSize
    int branchSize = (int)(atof(str.c_str()));
    
    while(ifs)
    {
        ifs >> str;
        CPnum = (int)(atof(str.c_str()));
        string firstFourBits = bitset<4>(CPnum).to_string();
       
        ifs >> str;
        degree = (int)(atof(str.c_str()));
        string lastFourBits = bitset<4>(degree).to_string();

        if(ifs.fail())  break;

        //of_uncompressed << firstFourBits << " - " << lastFourBits << '\n';
        of_uncompressed << CPnum << " - " << degree << '\n';
        writeBits(ofBuffer, firstFourBits + lastFourBits);

        ifs >> str;
        string SampleNum = bitset<16>((int)(atof(str.c_str()))).to_string();
        of_uncompressed << (int)(atof(str.c_str())) << '\n';
        writeBits(ofBuffer, SampleNum);

        while(CPnum--)  
        {
            ifs >> str;
            currentPoint[0] = (atof(str.c_str()));
            ifs >> str;
            currentPoint[1] = (atof(str.c_str()));
            ifs >> str;
            currentPoint[2] = (atof(str.c_str()));
            //ifs >> str;////4D
            //imp = (atof(str.c_str()));

            //////find---check if this point (its nerghbor) has already exit in the vector////
            if(storeCPs.empty())//first time
                storeCPs.push_back(currentPoint);
            else
            {
               for(auto it = storeCPs.begin();it!=storeCPs.end();it++)
               {
                   Vector3<int> point = *it;
                   if(abs(currentPoint[0]-point[0]) < 2 && abs(currentPoint[1]-point[1]) < 2 && abs(currentPoint[2]-point[2]) < 2)
                   {
                        RepeatNum ++;
                        find = true;
                        break;
                   }

               }
                if(!find)
                {
                    storeCPs.push_back(currentPoint);
                }
                find = false;  
            }
           ////end


            encode_x = currentPoint[0] - last_x;//maybe negive value.
            encode_y = currentPoint[1] - last_y;//maybe negive value.
            encode_dt = currentPoint[2] - last_dt;//maybe negive value.
            if (encode_x == 100 && encode_y == 100) {cout<<"repeat!!!!!"<<endl; encode_x = 101;}// To avoid conflict with the flag.
            if (encode_x > 127 || encode_y > 127 || encode_dt > 127 || encode_x < -127 || encode_y < -127 || encode_dt < -127)
            {
                //of_uncompressed << "0000000000000000" << encode_x << " / " << encode_y << " / " << encode_dt << '\n';
                //writeBits(ofBuffer, "0000000000000000");
                of_uncompressed << "0110010001100100" << encode_x << " / " << encode_y << " / " << encode_dt << '\n';
                writeBits(ofBuffer, "0110010001100100");
                
                writeBits(ofBuffer, encode_16bit_point(encode_x,encode_y,encode_dt));
            }
            else
            {
                of_uncompressed << encode_x << " / " << encode_y << " / " << encode_dt << '\n';
                writeBits(ofBuffer, encode_8bit_point(encode_x,encode_y,encode_dt));
            }
            last_x = currentPoint[0];
            last_y = currentPoint[1];
            last_dt = currentPoint[2];

        } 
    } 
    AddToOutput<<"The number of Repeat points is "<<RepeatNum<<endl;
    AddToOutput<<"The output size needs to be subtracted "<<branchSize*2<<" Bytes."<<endl;//since we don't need the SampleNum.
    AddToOutput.close();
    writeBits(ofBuffer, "11111111");
}

//original 16-bit method.

void ImageWriter::writeCP(unsigned int width, unsigned int height) {
    ifstream ifs("controlPoint.txt"); 
    string str;
    int x,y,dt,imp; ////if we use uint16_t, then we cannot represent negative value.(-5 turns to 65531)
    int CPnum, degree;
    
    num_pixels = width * height;

    of_uncompressed << width << " - " << height << '\n';
    writeBits(ofBuffer, bitset<16>(width).to_string() + bitset<16>(height).to_string());

    ifs >> str;//branchSize

    while(ifs)
    {
        ifs >> str;
        CPnum = (int)(atof(str.c_str()));
        string firstFourBits = bitset<4>(CPnum).to_string();
       
        ifs >> str;
        degree = (int)(atof(str.c_str()));
        string lastFourBits = bitset<4>(degree).to_string();

        if(ifs.fail())  break;

        //of_uncompressed << firstFourBits << " - " << lastFourBits << '\n';
        of_uncompressed << CPnum << " - " << degree << '\n';
        writeBits(ofBuffer, firstFourBits + lastFourBits);

        ifs >> str;
        string SampleNum = bitset<16>((int)(atof(str.c_str()))).to_string();
        of_uncompressed << (int)(atof(str.c_str())) << '\n';
        writeBits(ofBuffer, SampleNum);

        while(CPnum--)  
        {
            ifs >> str;
            x = (atof(str.c_str()));
            ifs >> str;
            y = (atof(str.c_str()));
            ifs >> str;
            dt = (atof(str.c_str()));
            //ifs >> str;////4D
            //imp = (atof(str.c_str()));

            of_uncompressed << x << " - " << y << " - " << dt << '\n';
            
            writeBits(ofBuffer, encode_16bit_point(x,y,dt));
        } 
    }    
    writeBits(ofBuffer, "11111111");
}


