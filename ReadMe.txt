This software package is the technical support for our paper "Spline-based medial axis transform representation of binary images". 

# 1. Building


The software needs a C++ compiler, OpenGL, the CUDA, and GLUT to build. 
Under a Debian based linux they can be installed using the following command.

* install

sudo apt-get install cuda zpaq libzstd-dev libsnappy-dev gcc git make cmake libboost-all-dev freeglut3-dev libglew-dev ragel libvala-0.40-dev glmark2 valac liblz4-dev liblzma-dev libbz2-1.0

* to build:

cd Code/imShow && make

cd ../../

cd Code/imConvert && make


# 2. Running


cd Code/imConvert

./skeletonify config.txt

cd ../imShow

./show_skeleton output.sir 

where config.txt is a configuration file. The detailed explanation is shown in the last section.


**I also write a bash script to execute the whole pipeline, named 'pipeline' in the Code file.**
What you need to do are:
 
(a) Put binary images that you want to test into imConvert/DATA. 

(b) Change parameters as you want in this script. These parameters include what saliency threshold do you need ("Saliency"), what's the Hausdorff distance between skeleton branches and the corresponding fitted splines do you need ("hausdorff"), Do you need merge procedure or not ("merge" is set to 1 or 0), etc, see the detailed explanation in the 4th section.

(c) Run ./pipeline. Then all quality scores and various sizes data will be written into 'output.txt'. Besides, in computeError file, you'll see the generated result - output.png, in which cyan shows the SMAT-based reconstruction of the shape while the black outline represents the original shape. The fitted splines and their corresponding control points are also shown in this image.


# 3. Function


imConvert
---------

This program reads a binary image (a PGM image in the DATA file), extracts its skeleton, segments the skeleton to several branches and then fitted with splines. Next, the control points of the fitted spline will be saved to the controlPoint.txt. Then im.writeCP will read these control points and then saved in into the 'output.smat' file which located in the imShow file. 

SMAT format is:

UseZPAQorNot (8bit) - Always set to 0.

WIDTH (16bit) - The width of the image.

HEIGHT (16bit) - The height of the image.

BranchCPdata - Control points of branches are stored here. Size varies.

END (8bit) - Use 11111111 as the end tag.

For each branch, BranchCPdata format is like this:

CPnum (4bit) + degree (4bit) - Both degree and the number of CP are necessary for the generation of a spline, and both of them are less than 16. So to save space, I use one byte to store them. The high nibble represents the number of control points of this branch, and it determines how many bytes should be read next for this branch.

SampleNum (16bit) - The sample number of this branch.

CPdata (6*CPnum bytes) - CPdata is stored as: X (16bit) - The initial X coordinate. Y (16bit) - The initial Y coordinate.  DT (16bit) - The initial DT value.


imShow
---------
This program is responsible for decode the SMAT file and then do the reconstruction.

* getAlphaMapOfLayer() - This function draws all disks. The alpha will be 1 at the location that should be drawn, otherwise 0. 


# 4. Other remarks

As said before, config.txt is a configuration file. The following explains all parameters:

* outputFile 

This parameter indicates where to put the output file, i.e., the .smat file.

* SkeletonThreshold 

This parameter is a support for the saliency skeleton. After the saliency skeleton thresholding, there will be a collection of disconnected, skeletal components. Only the longest 'core' component/fragment should be preserved. So this parameter is used to filter those small disconnected skeleton branches. See details in the paper "Feature Preserving Smoothing of Shapes using Saliency Skeletons".

* Imp

Threshold for Importance metric. See details in the paper "Feature Preserving Smoothing of Shapes using Saliency Skeletons". Set it to 3 is a good one. 

* removeRepCR

It is generally set to 0, means do not remove repeated point, and use 8-bit encoding scheme.

* filename

Input image, must be PGM file.

* ssThreshold

An important parameter. It represents saliency skeleton threshold. 0.2-2 should be a good range.

* hausdorff

Another important parameter. It tells how accurately B-splines fit medial branches. The smaller this error, the higher accuracy. 0.001 - 0.005 should be a good range.

