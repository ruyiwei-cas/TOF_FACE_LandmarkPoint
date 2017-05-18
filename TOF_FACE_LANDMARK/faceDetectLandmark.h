#pragma once
#include <cstdio>
#include <opencv2/opencv.hpp>
#include <fstream>

using namespace std;
using namespace cv;

// dirty but works
//int train(int);
//int test(void);
//int prepare(void);
int run(Mat img, vector<CvRect>& rect, vector<Mat>& landmarkPoint, const string modelPath);
