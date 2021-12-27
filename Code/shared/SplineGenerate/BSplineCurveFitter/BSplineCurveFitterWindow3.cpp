// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#include "BSplineCurveFitterWindow3.h"
#include <Graphics/VertexColorEffect.h>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
  
 
using namespace std;
#define mDimension 3
#define colorNum 10000

int color[colorNum][3];//wang
 
int num = 0, CPNum = 0, cntlSepNum=0,  orignum=0;//how many parts you should draw
std::shared_ptr<Visual> OrigCurve[10000];//if we put this in .h file, it would cause the error of "corrupted size vs. prev_size". I donot know why.
std::shared_ptr<Visual> drawSpline[10000];//give enough space first
std::shared_ptr<Visual> drawControlPoint[1000];
std::shared_ptr<Visual> drawControlPointl[10000];
float maxhaus = 0.0f;
int maxBranch = 0;
int totalSample = 0;
int init = 1;//which branch seq you should start with.
int TotalControlNum = 0; 
int DeterminedNumControls, DeterminedDegree, lastNumCntl;
float minError; bool connect=0;
vector<Vector3<float>> controlData(10000);
int countSepCntl = 0;
float lastx, lasty;
ofstream OutFile, OutFile1, OutControlNum, outMerge, OutFileBranch;
int width, height;
float minErrorThreshold;
bool gapfilling;
//vector<vector<Vector3<float>>> sampleSet(1000);
float storevector[mDimension]={0,0,0};
vector<Vector3<float>> controlDataVec(1000);
unique_ptr<BSplineCurveGenerate<float>> SplineGeneratePtr;
int mergeOrNot = true;
bool deleteshort = false;
int splitNum = 0;
unsigned int MinAllowableLength = 5;
int branchNum = 0;
/*//
BSplineCurveFitterWindow3::BSplineCurveFitterWindow3(vector<vector<Vector3<float>>> BranchSet, float hausdorff_,float diagonal_)
    :
    sampleSet(BranchSet)
{
    minErrorThreshold = hausdorff_;
    diagonal = diagonal_;
    cout<<"good-----!!"<<endl;
    SplineFit();
    //generateColor();
    //CreateScene();
    //InitializeCamera(60.0f, GetAspectRatio(), 0.1f, 100.0f, 0.01f, 0.001f,
    //    { 0.0f, 0.0f, -4.0f }, { 0.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 0.0f });
    //mPVWMatrices.Update();
}
*/
void readConfig()
{
    string str1;
    ifstream ifs("/home/jieying/Desktop/SMAT/Code/imShow/config.txt");
    ifs >> str1;
    mergeOrNot = (int)atof(str1.c_str());
}

BSplineCurveFitterWindow3::BSplineCurveFitterWindow3()
{
    //cout<<"good-----!!"<<endl;
    readConfig();
}

int BSplineCurveFitterWindow3::SplineFit(vector<vector<Vector3<float>>> BranchSet, float hausdorff_,float diagonal_, vector<int *> connection_)
{
    sampleSet = BranchSet;
    minErrorThreshold = hausdorff_;
    diagonal = diagonal_;
    connection = connection_;
//cout<<"-----------------"<<sampleSet.size()<<endl;
    if (mergeOrNot) {NewMerge();cout<<"new"<<endl;}
    //else if (mergeOrNot==2) Merge();
//cout<<"-----------------"<<sampleSet.size()<<endl;
    OutFile1.open("controlPoint.txt");
    OutControlNum.open("../output.txt",ios_base::app);
    int sampleSetSizeMinus = 0;

    for (unsigned int i = 0; i < sampleSet.size();i++)
    {
        if (sampleSet[i].size()<MinAllowableLength) 
            sampleSetSizeMinus++;
    }
    OutFile1 << (sampleSet.size()-sampleSetSizeMinus)<<endl;
clock_t start = clock();
    for (unsigned int i = 0; i < sampleSet.size();i++)
    {
        if (sampleSet[i].size()<MinAllowableLength) continue;
        CreateBSplinePolyline(sampleSet[i]);
    }    
    clock_t endtime = clock();
    cout<<"run time: "<<(double)(endtime-start)/CLOCKS_PER_SEC<<" s" <<endl;

    OutControlNum<<"TotalControlNum: "<<TotalControlNum<<endl;
    OutFile1<<"30000"<<endl; //just a sign
    OutFile1.close();  
    
    OutControlNum.close();
    return splitNum;
}

void BSplineCurveFitterWindow3::SplineGenerate(int SuperR)
{
    int CPnum,degree;
    unsigned int numSamples;
    string str;
    float controlData;
    OutFile.open("sample.txt");
    OutFileBranch.open("BranchSample.txt");

    vector<float> mControlData;

   ifstream ifs1("ControlPoints.txt"); 
    //ifstream ifs1("1.txt"); 
    ifs1 >> str;
    int width = (int)atof(str.c_str());
    OutFile << width <<" ";
    ifs1 >> str;
    int height = (int)atof(str.c_str());
    OutFile << height <<endl;

    diagonal = sqrt((float)(width * width + height * height));
    
    while(ifs1)
    { 
        //count1++;
        
        ifs1 >> str;
        CPnum = (int)atof(str.c_str());
        ifs1 >> str;
        degree = (int)atof(str.c_str());
        ifs1 >> str;
        numSamples = (unsigned int)atof(str.c_str());
    
        if(ifs1.fail())  break;
        
        for (int i = 0; i< CPnum; ++i)
        {
            for (int j = 0; j < mDimension; ++j)
            {
                ifs1 >> str;
                controlData = atof(str.c_str())/diagonal;
                controlDataVec[i][j] = controlData;
                mControlData.push_back(controlData); 
            }
        }
            
            
        TotalControlNum += CPnum;
        totalSample += numSamples;
        SplineGeneratePtr = std::make_unique<BSplineCurveGenerate<float>>(mDimension, degree, mControlData, CPnum);
        //cout<<"mControlData.size(): "<<mControlData.size()/3<<" count: "<<count1<<endl;
       
        controlDataVec.resize(CPnum);
        CreateGraphics(numSamples, SuperR);
        mControlData.clear();
    }
    OutFileBranch.close();
    OutFile.close(); 
    
    ifs1.close();

}


void BSplineCurveFitterWindow3::CreateGraphics(unsigned int numSamples, int superR)
{
    branchNum++;
    unsigned int numSplineSample = (unsigned int)(numSamples*2*superR);//sub-pixel.
   // unsigned int numSplineSample = numSamples; //uniform sampling
    float vector[mDimension]={0,0,0};
    float multiplier = 1.0f / (numSplineSample - 1.0f);

    for (unsigned int i = 0; i < numSplineSample; ++i)
    {
        float t = multiplier * i;
       
        SplineGeneratePtr->GetPosition(t, reinterpret_cast<float*>(vector));
        OutFile<<vector[0]*diagonal<<" "<<vector[1]*diagonal<<" "<<vector[2]*diagonal<<endl;      //save to the txt file.
        OutFileBranch<<branchNum<<" "<<vector[0]*diagonal<<" "<<vector[1]*diagonal<<endl;
    }
}

float BSplineCurveFitterWindow3::Judge(vector<Vector3<float>> Sample)
{
    unsigned int numSamples = (unsigned int)Sample.size();
    vector<Vector3<float>> SplineSamples(10000);
    
    float CPandError = 0;
    
    unsigned int numSplineSamples = (unsigned int)(numSamples * 1.4);
    //unsigned int numSplineSamples = numSamples; // uniform sampling.
    float multiplier = 1.0f / (numSplineSamples - 1.0f);
    SplineSamples.resize(numSplineSamples);
    
    int minDegree;
    int numControls = 2;
   // int iter = 0;
    while (1)
    {
        //iter ++;
        minError = 100.0f; minDegree = 10;
        for (int degree = 1; degree < numControls; degree++)
        {
            mSpline = std::make_unique<BSplineCurveFit<float>>(mDimension, static_cast<int>(Sample.size()),
            reinterpret_cast<float const*>(&Sample[0]), degree, numControls);

            for (unsigned int i = 0; i < numSplineSamples; ++i)
            {
                float t = multiplier * i;
                mSpline->GetPosition(t, reinterpret_cast<float*>(storevector));
                 for(int y=0;y<mDimension;y++)
                    SplineSamples[i][y] = storevector[y];
            }
            
        // Compute error measurements.
            float maxLength = 0.0f;
            Vector3<float> diff;
            float sqrLength, minLength;
            for (unsigned int i = 0; i < numSamples; ++i)
            {
                minLength = 100.0f;
                for (unsigned int j = 0; j < numSplineSamples; ++j)
                {
                    diff = Sample[i] - SplineSamples[j];
                    sqrLength = Dot(diff, diff);
                    if (sqrLength < minLength) minLength = sqrLength;
                }
                if (minLength > maxLength) maxLength = minLength;
            
            }
            hausdorff = std::sqrt(maxLength);
            if (minError > hausdorff) { minError = hausdorff; minDegree = degree;}
            //cout<<numControls<<" hausdorff: "<<hausdorff<<" degree: "<<degree<<endl;
        }
        //cout<<numControls<<" minError: "<<minError<<" minDegree: "<<minDegree<<endl;
        if (numControls > 12) //means this branch is not easy to fit well, so I split it two half.
        {
            CPandError = 100;//assign an big enough value
            return CPandError;
            //break;
        }
        if (minError < minErrorThreshold)
        {
            DeterminedNumControls = numControls;
            DeterminedDegree = minDegree;
            //cout<<" minError: "<<minError<<endl;
            CPandError = numControls + minError;

            break;
        }
        else numControls ++;
    }
   
    return CPandError;
}

void BSplineCurveFitterWindow3::CreateBSplinePolyline(vector<Vector3<float>> Sample)
{
    float cpError = Judge(Sample);
    
    if (cpError == (float)100) //the branch may be too long to fit well.
    {
        cout<<"---------"<<endl;
        splitNum ++;
        vector<Vector3<float>> first = Sample;//init to VOID segmentation fault.
        vector<Vector3<float>> second = Sample;
        
        for (unsigned int i = Sample.size()/2; i < Sample.size(); ++i)
            first[i] = 0;
        first.resize(Sample.size()/2);  
        for (unsigned int i = Sample.size()/2; i < Sample.size(); ++i)
            second[i - Sample.size()/2] = Sample[i];
        second.resize(Sample.size() - Sample.size()/2);  
   
        CreateBSplinePolyline(first);
        CreateBSplinePolyline(second);
    }
    else{
        mSpline = std::make_unique<BSplineCurveFit<float>>(mDimension, static_cast<int>(Sample.size()),
            reinterpret_cast<float const*>(&Sample[0]), DeterminedDegree, DeterminedNumControls);
   
        TotalControlNum += DeterminedNumControls;

        OutFile1<< DeterminedNumControls <<" "<<DeterminedDegree<<" "; 
        OutFile1<< Sample.size() <<" ";  

        float const* controlDataPtr = mSpline->GetControlData();
        for (int i = 0; i< DeterminedNumControls; ++i)
        {
            for (int j = 0; j < mDimension; ++j)
            {
                controlData[i][j] = (*controlDataPtr);
                OutFile1<<(round)(*controlDataPtr*diagonal)<<" ";
                //mControlData.push_back(*controlDataPtr);
                controlDataPtr++;
            }
            
        }    
        OutFile1<<endl;
    }  
}

void BSplineCurveFitterWindow3::NewMerge()
{
    vector<Vector3<float>> first, second;
    float minEandContlNum = 100.0;
    //float maxdiff = 0.0;
    int minIndex = 1000;
    float firstE,secondE,mergeE;
    
    //vector<int> GapFill;////
    //outMerge.open("outMerge.txt");
  //time:0.007
    //clock_t start = clock();
    for(auto it = connection.begin();it!=connection.end();it++)
    {
        int *sampleIndex = *it;
        //cout<<out[0]<<"--"<<out[1]<<"--"<<out[2]<<"--"<<out[3]<<endl;
        first = sampleSet[sampleIndex[0]];
        if (deleteshort) firstE = Judge(first);
        else 
        {
            if (first.size()<MinAllowableLength) firstE = 2.01;//
            else firstE = Judge(first);
        }

        for(int index = 1; index < 4; index++)//index for 'sampleIndex'.
        {
            if(sampleIndex[index]==0) continue;
            second = sampleSet[sampleIndex[index]];
                
            if (deleteshort) secondE = Judge(second);
            else 
            {
                if(second.size()<MinAllowableLength) secondE = 2.01;//assign a big num.
                else secondE = Judge(second);
            }

            merge.resize(10000);//This is very important!!
            for (unsigned int i = 0; i < first.size(); ++i)
                merge[i] = first[i];
            
            for (unsigned int i = first.size(); i < (first.size()+second.size()); ++i)
                merge[i] = second[i-first.size()];
                
            merge.resize(first.size()+second.size());
           

            if (deleteshort) mergeE = Judge(merge);
            else
            {
                if (merge.size()<MinAllowableLength) mergeE = 4; //about 3CP + 0.- error.
                else mergeE = Judge(merge);
            }
    
            //cout<<"mergeE: "<<mergeE<<" firstE: "<<firstE<<" secondE: "<<secondE<<endl;
            if (mergeE < (firstE + secondE)) 
                if (minEandContlNum > mergeE) { minEandContlNum = mergeE; minIndex = sampleIndex[index];}
            
        }
        if (minIndex!=1000) //meet the merge condition.
        {
            //cout<<"Merge！ first: "<<sampleIndex[0]<<" second: "<<minIndex<<endl;
            second = sampleSet[minIndex];
            merge.resize(10000);//This is very important!!
            for (unsigned int i = 0; i < first.size(); ++i)
                merge[i] = first[i];
            
            for (unsigned int i = first.size(); i < (first.size()+second.size()); ++i)
                merge[i] = second[i-first.size()];
            merge.resize(first.size()+second.size());

            //for (unsigned int i = 0; i<merge.size(); i++)
            //    outMerge<<iter<<" "<<merge[i][0]*diagonal<<" "<<merge[i][1]*diagonal<<" "<<merge[i][2]*diagonal<<endl;
        

            sampleSet.erase(sampleSet.begin()+sampleIndex[0]);
            sampleSet.insert(sampleSet.begin()+sampleIndex[0],merge);
            //sampleSet.erase(sampleSet.begin()+minIndex);
            sampleSet[minIndex].clear();//still occupy the position.

            minEandContlNum = 100.0f;
            minIndex = 1000;
        }
        
    }
    
    cout<<"--"<<sampleSet.size()<<endl;
    for (unsigned int i = sampleSet.size()-1; i > 0;i--)
    {
        if (sampleSet[i].size()<1) 
            sampleSet.erase(sampleSet.begin()+i);
    }


    
    //outMerge.close();
    cout<<"--"<<sampleSet.size()<<endl;
     
}

//original, slow algorithm

void BSplineCurveFitterWindow3::Merge()
{
    vector<Vector3<float>> first, second;
    float minEandContlNum = 100.0;
    //float maxdiff = 0.0;
    int minIndex;
    float firstE,secondE,mergeE;
    int iter = 0;
    vector<int> GapFill;////
    outMerge.open("outMerge.txt");
  //time:0.007
    //clock_t start = clock();
    bool finish = 0;
    while (!finish)
    {
        iter++;   
        for (unsigned int i = 0; i<sampleSet.size();i++)
        {
            GapFill.clear();//important
            first = sampleSet[i];
            if (deleteshort) firstE = Judge(first);
            else 
            {
                if (first.size()<MinAllowableLength) firstE = 2.01;//
                else firstE = Judge(first);
            }
            for (unsigned int j = i+1;j<sampleSet.size();j++)
            {
                second = sampleSet[j];
                //cout<<"i: "<<i<<" j: "<<j<<" - "<<second[0][0]*diagonal<<" - "<<second[0][1]*diagonal<<endl;
                if (abs(first[first.size()-1][0]-second[0][0]) < 2.0/diagonal && abs(first[first.size()-1][1]-second[0][1]) < 2.0/diagonal)//neighbor  //add other situations
                {
                    //cout<<"nerghbor: "<<iter<<" i: "<<i<<" j: "<<j<<endl;
                    GapFill.push_back(j);

                    if (deleteshort) secondE = Judge(second);
                    else 
                    {
                        if(second.size()<MinAllowableLength) secondE = 2.01;//assign a big num.
                        else secondE = Judge(second);
                    }

                    merge.resize(10000);//This is very important!!
                    for (unsigned int i = 0; i < first.size(); ++i)
                        merge[i] = first[i];
                    /*    //remove repeat junction points
                    if (abs(first[first.size()-1][0]-second[0][0]) < 1.0/diagonal && abs(first[first.size()-1][1]-second[0][1]) < 1.0/diagonal)//repeat points 
                    {
                        for (unsigned int i = first.size(); i < (first.size()+second.size())-1; ++i)
                            merge[i] = second[i-first.size()+1];
                        merge.resize(first.size()+second.size()-1);//important
                    }
                    else
                    {
                        for (unsigned int i = first.size(); i < (first.size()+second.size()); ++i)
                            merge[i] = second[i-first.size()];
                        merge.resize(first.size()+second.size());//important
                
                    }
                    */
                   for (unsigned int i = first.size(); i < (first.size()+second.size()); ++i)
                        merge[i] = second[i-first.size()];
                        
                    merge.resize(first.size()+second.size());
                    //for (unsigned int i = 0; i<merge.size(); i++)
                    //    outMerge<<iter<<" "<<merge[i][0]*diagonal<<" "<<merge[i][1]*diagonal<<" "<<merge[i][2]*diagonal<<endl;


                    if (deleteshort) mergeE = Judge(merge);
                    else
                    {
                        if (merge.size()<MinAllowableLength) mergeE = 4; //about 3CP + 0.- error.
                        else mergeE = Judge(merge);
                    }
            
                    //cout<<"mergeE: "<<mergeE<<" firstE: "<<firstE<<" secondE: "<<secondE<<endl;
                    if (mergeE < (firstE + secondE)) 
                        if (minEandContlNum > mergeE) { minEandContlNum = mergeE; minIndex = j;}
                    
                }
            }
            if (minIndex!=100) //meet the merge condition.
            {
                //cout<<"Merge！ i: "<<i<<" j: "<<minIndex<<endl;
                second = sampleSet[minIndex];
                merge.resize(10000);//This is very important!!
                for (unsigned int i = 0; i < first.size(); ++i)
                    merge[i] = first[i];
                /* //remove repeat junction points
                if (abs(first[first.size()-1][0]-second[0][0]) < 1.0/diagonal && abs(first[first.size()-1][1]-second[0][1]) < 1.0/diagonal)//repeat points 
                {
                    for (unsigned int i = first.size(); i < (first.size()+second.size())-1; ++i)
                        merge[i] = second[i-first.size()+1];
                    merge.resize(first.size()+second.size()-1);
                }
                else
                {
                    for (unsigned int i = first.size(); i < (first.size()+second.size()); ++i)
                        merge[i] = second[i-first.size()];
                    merge.resize(first.size()+second.size());
                }
                */
                for (unsigned int i = first.size(); i < (first.size()+second.size()); ++i)
                    merge[i] = second[i-first.size()];
                merge.resize(first.size()+second.size());

                for (unsigned int i = 0; i<merge.size(); i++)
                    outMerge<<iter<<" "<<merge[i][0]*diagonal<<" "<<merge[i][1]*diagonal<<" "<<merge[i][2]*diagonal<<endl;
            

                sampleSet.erase(sampleSet.begin()+i);
                sampleSet.insert(sampleSet.begin()+i,merge);
                sampleSet.erase(sampleSet.begin()+minIndex);

                minEandContlNum = 100.0f;
                minIndex = 100;
                break;
            }
            else if (i==sampleSet.size()-1) finish = 1;
        }

    } //finished.
    //clock_t endtime = clock();
//cout<<"run time: "<<(double)(endtime-start)/CLOCKS_PER_SEC<<" s" <<endl;
    outMerge.close();
    cout<<"--"<<sampleSet.size()<<endl;
     
}