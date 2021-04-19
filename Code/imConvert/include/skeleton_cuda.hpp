#ifndef SKEL_CUDA_HPP
#define SKEL_CUDA_HPP

FIELD<float>* computeSkeleton(FIELD<float> *im);
int initialize_skeletonization(FIELD<float>* im);
void analyze_cca(FIELD<float>* skel);
void deallocateCudaMem();
short* get_current_skel_ft();
FIELD<float>* skelft_to_field();

#endif