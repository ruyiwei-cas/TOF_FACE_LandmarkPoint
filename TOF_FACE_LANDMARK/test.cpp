#include "./lbf/lbf.hpp"
#include "faceDetect.h"
#include <cstdio>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
using namespace std;
using namespace lbf;

void parseTxt(string &txt, vector<Mat> &imgs, vector<Mat> &gt_shapes, vector<BBox> &bboxes);

int run(Mat img,vector<CvRect>& rects,vector<Mat>& landmarkPoint,const string modelPath) {
    Config &config = Config::GetInstance();
    int landmark_n = config.landmark_n;
 
    LbfCascador lbf_cascador;
    FILE *model = fopen(config.saved_file_name.c_str(), "rb");
    lbf_cascador.Read(model);
    fclose(model);

	CvSize cvSize; cvSize.width = img.cols; cvSize.height = img.rows;
	IplImage *imgIplImage = cvCreateImage(cvSize, 8, img.channels()); 
	memcpy(imgIplImage->imageData, (char*)img.data, sizeof(char)*img.channels()*img.rows*img.cols);
	detectMulti(imgIplImage, rects, modelPath);
	
	for (int i = 0; i < rects.size(); i++)
	{
		double x_, y_, w_, h_;
		x_ = 0; y_ = img.rows/2;
		w_ = img.cols; h_ = img.rows/2;
		BBox bbox_(rects[i].x - x_, rects[i].y - y_, rects[i].width,rects[i].height);
		Rect roi(x_, y_, w_, h_);
		img = img(roi).clone();

		Mat gray;
		if (img.channels() == 3)
			cvtColor(img, gray, CV_BGR2GRAY);
		else
			gray = img;
		Mat shape = lbf_cascador.Predict(gray, bbox_);
		landmarkPoint.push_back(shape);
#ifdef _DEBUG
		img = drawShapeInImage(img, shape, bbox_);
		resize(img, img, Size(640, 640));
		imshow("landmark", img);
		waitKey(0);
#endif 
	}
  
    return 0;
}
