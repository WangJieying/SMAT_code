/*
 * File:   ImageWriter.hpp
 * Author: yuri
 *
 * Created on June 6, 2011, 2:59 PM
 */

#ifndef IMAGEWRITER_HPP
#define IMAGEWRITER_HPP

#include <sstream>
#include <vector>
#include "CUDASkel2D/include/field.h"
#include "Triple.hpp"
#include "SkelNode.hpp"
#include <math.h>
// #include <boost/dynamic_bitset.hpp>
// #include "../shared/FastAC/arithmetic_codec.h"
#include "ImageEncoder.hpp"
// #include <utility>
// #include <numeric>

using namespace std;

typedef std::pair<int, int> coord2D_t;
typedef vector<coord2D_t> coord2D_list_t; 
typedef Triple<int, int, int> coord3D_t;
typedef SkelNode< coord3D_t > skel_tree_t;
typedef boost::dynamic_bitset<uint8_t> BitSet;
typedef pair<int, pair<int, int>> loc_ding;

class ImageWriter {
public:
    ImageWriter(const char *fname);
    void save();
    virtual ~ImageWriter();
   
    //void writeHeader(unsigned int width, unsigned int height);
    void writeCP(unsigned int width, unsigned int height);
    void writeCP_(unsigned int width, unsigned int height);
    void writeCPandID(unsigned int width, unsigned int height);
    
    //void write_color_image(vector<std::pair<int, skel_tree_t*>>* red_forest, vector<std::pair<int, skel_tree_t*>>* green_forest, vector<std::pair<int, skel_tree_t*>>* blue_forest);

private:
    /*******************************************************************************/
    /*******************************************************************************/
    /** IMPORTANT -- WHEN CHANGES ARE MADE TO THE FILE FORMAT, UPDATE THIS NUMBER **/
    uint16_t WRITER_FILE_VERSION_NUMBER = 0xB;
    uint16_t READER_FILE_VERSION_NUMBER = 0xB;
    /*******************************************************************************/
    /*******************************************************************************/
    enum COMPRESS_MODE {COMPRESS_LZMA = 0, COMPRESS_ZLIB, COMPRESS_LZHAM, COMPRESS_LZMA2, COMPRESS_BROTLI, COMPRESS_ZPAQ, COMPRESS_BSC, COMPRESS_CSC, COMPRESS_BZIP2, COMPRESS_ZSTD};
    enum ENCODE_MODE {HUFFMAN = 1, CANONICAL, UNITARY, EXP_GOULOMB, ARITHMETIC, COMPACT, MTF, PREDICTIVE, RAW, TRADITIONAL};
   
    void writeBits(ostream& os, const string& str);
   
    ImageEncoder* encoder;

    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> mark_important_points(vector<std::pair<int, skel_tree_t*>>* forest);
    map<pair<int, int>, vector<pair<int, int>>> flatten(vector<std::pair<int, skel_tree_t*>>* forest, bool print = false);
    // void bundle_points(vector<std::pair<int, skel_tree_t*>>* forest);
    uint8_t out8[3];
    uint16_t width, height;
    uint num_pixels;
    stringstream ofBuffer;
    ofstream of;
    ofstream of_uncompressed;
    COMPRESS_MODE comp_mode;
    ENCODE_MODE encode_mode;
    std::string codec_name;
    const char *fname;
    typedef struct {
        COMPRESS_MODE mode;
        string compress_algo;
        string in_data;             // The data that is to be compressed
        size_t in_size;             // size of the uncompressed data
        size_t dest_len;            // compressed size
        unsigned char* out_data;    // compressed data

    } ret_struct;
    string get_compression_name(COMPRESS_MODE mode);

    void set_compression_algorithm(std::string codec);

    /* Statistics */
    int NUMPOINTS, NUMPATHS, OLDNUMPOINTS, OLDNUMPATHS;
  
};


#endif  /* IMAGEWRITER_HPP */

