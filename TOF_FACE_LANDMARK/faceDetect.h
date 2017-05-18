#pragma once
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include "picornt.h"
typedef unsigned char  uint8_t;



void process_image(IplImage* frame, std::vector<CvRect>& rect, int draw);
int detectMulti(IplImage* img, std::vector<CvRect>& faceRect, std::string modelPath);
