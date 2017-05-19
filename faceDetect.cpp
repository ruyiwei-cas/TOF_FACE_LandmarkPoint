#include"faceDetect.h"
#ifdef __GNUC__
#include <time.h>
float getticks()
{
	struct timespec ts;

	if (clock_gettime(CLOCK_MONOTONIC, &ts) < 0)
		return -1.0f;

	return ts.tv_sec + 1e-9f*ts.tv_nsec;
}
#else
#include <windows.h>

float getticks()
{
	static double freq = -1.0;
	LARGE_INTEGER lint;

	if (freq < 0.0)
	{
		if (!QueryPerformanceFrequency(&lint))
			return -1.0f;

		freq = lint.QuadPart;
	}

	if (!QueryPerformanceCounter(&lint))
		return -1.0f;

	return (float)(lint.QuadPart / freq);
}
#endif


void* cascade = 0;

int minsize;
int maxsize;

float angle;

float scalefactor;
float stridefactor;

float qthreshold;

int usepyr;
int noclustering;
int verbose;
#define MAXNDETECTIONS 5

void process_image(IplImage* frame,std::vector<CvRect>& rect, int draw)
{
	int i, j;
	float t;

	uint8_t* pixels;
	int nrows, ncols, ldim;

	int ndetections;
	float qs[MAXNDETECTIONS], rs[MAXNDETECTIONS], cs[MAXNDETECTIONS], ss[MAXNDETECTIONS];

	IplImage* gray = 0;
	IplImage* pyr[5] = { 0, 0, 0, 0, 0 };
	
	gray = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, 1);
	if (!pyr[0]&& usepyr)
	{//对图像进行金字塔的缩放
		pyr[0] = gray;
		pyr[1] = cvCreateImage(cvSize(frame->width / 2, frame->height / 2), frame->depth, 1);
		pyr[2] = cvCreateImage(cvSize(frame->width / 4, frame->height / 4), frame->depth, 1);
		pyr[3] = cvCreateImage(cvSize(frame->width / 8, frame->height / 8), frame->depth, 1);
		pyr[4] = cvCreateImage(cvSize(frame->width / 16, frame->height / 16), frame->depth, 1);
	}

	// get grayscale image
	if (frame->nChannels == 3)
		cvCvtColor(frame, gray, CV_RGB2GRAY);
	else
		cvCopy(frame, gray, 0);

	// perform detection with the pico library
	t = getticks();

	if (usepyr)
	{
		int nd;

		//
		pyr[0] = gray;

		pixels = (uint8_t*)pyr[0]->imageData;
		nrows = pyr[0]->height;
		ncols = pyr[0]->width;
		ldim = pyr[0]->widthStep;

		ndetections = find_objects(rs, cs, ss, qs, MAXNDETECTIONS, cascade, angle, pixels, nrows, ncols, ldim, scalefactor, stridefactor, MAX(16, minsize), MIN(128, maxsize));

		for (i = 1; i < 5; ++i)
		{
			cvResize(pyr[i - 1], pyr[i], CV_INTER_LINEAR);

			pixels = (uint8_t*)pyr[i]->imageData;
			nrows = pyr[i]->height;
			ncols = pyr[i]->width;
			ldim = pyr[i]->widthStep;

			nd = find_objects(&rs[ndetections], &cs[ndetections], &ss[ndetections], &qs[ndetections], MAXNDETECTIONS - ndetections, cascade, angle, pixels, nrows, ncols, ldim, scalefactor, stridefactor, MAX(64, minsize >> i), MIN(128, maxsize >> i));

			for (j = ndetections; j < ndetections + nd; ++j)
			{//重新回调到具体的倍数
				rs[j] = (1 << i)*rs[j];
				cs[j] = (1 << i)*cs[j];
				ss[j] = (1 << i)*ss[j];
			}

			ndetections = ndetections + nd;
		}
	}
	else
	{
		//
		pixels = (uint8_t*)gray->imageData;
		nrows = gray->height;
		ncols = gray->width;
		ldim = gray->widthStep;
		//cs,rs 中心点的位置，ss大小，qs判断某个是否是图片的阈值，cascade模型，angle是寻找角度，pixels是指图像本身
		ndetections = find_objects(rs, cs, ss, qs, MAXNDETECTIONS, cascade, angle, pixels, nrows, ncols, ldim, scalefactor, stridefactor, minsize, MIN(nrows, ncols));
	}

	if (!noclustering)
		ndetections = cluster_detections(rs, cs, ss, qs, ndetections);

	CvRect faceBox;
	for (i = 0; i < ndetections; ++i){
		/*std::cout << "qs[i] = "<<qs[i] << std::endl;*/
		if (qs[i] >= qthreshold) {// check the confidence threshold
			faceBox.x = cs[i] - ss[i] / 2;
			faceBox.y = rs[i] - ss[i] / 2;
			faceBox.width = ss[i];
			faceBox.height = ss[i];
			rect.push_back(faceBox);
			cvRectangle(frame, cvPoint(cs[i] - ss[i] / 2, rs[i] - ss[i] / 2), cvPoint(cs[i] + ss[i] / 2, rs[i] + ss[i] / 2), CV_RGB(255, 0, 0), 4, 8, 0);
		/*	cvShowImage("faceImage", frame);
			cvWaitKey(0);*/
		}
	}

	for (int i = 0;i<5;i++){
		cvReleaseImage(&pyr[i]);
	}

}

//void process_webcam_frames()
//{
//	CvCapture* capture;
//
//	IplImage* frame;
//	IplImage* framecopy;
//
//	int stop;
//
//	const char* windowname = "--------------------";
//
//	//
//	capture = cvCaptureFromCAM(0);
//
//	if (!capture)
//	{
//		printf("* cannot initialize video capture ...\n");
//		return;
//	}
//
//	// the main loop
//	framecopy = 0;
//	stop = 0;
//
//	while (!stop)
//	{
//		// wait 5 miliseconds
//		int key = cvWaitKey(5);
//
//		// get the frame from webcam
//		if (!cvGrabFrame(capture))
//		{
//			stop = 1;
//			frame = 0;
//		}
//		else
//			frame = cvRetrieveFrame(capture, 1);
//
//		// we terminate the loop if the user has pressed 'q'
//		if (!frame || key == 'q')
//			stop = 1;
//		else
//		{
//			// we mustn't tamper with internal OpenCV buffers
//			if (!framecopy)
//				framecopy = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, frame->nChannels);
//			cvCopy(frame, framecopy, 0);
//
//			// webcam outputs mirrored frames (at least on my machines)
//			// you can safely comment out this line if you find it unnecessary
//			cvFlip(framecopy, framecopy, 1);
//
//			// ...
//			process_image(framecopy, 1);
//
//			// ...
//			cvShowImage(windowname, framecopy);
//		}
//	}
//
//	// cleanup
//	cvReleaseImage(&framecopy);
//	cvReleaseCapture(&capture);
//	cvDestroyWindow(windowname);
//}



int detectMulti(IplImage* img,std::vector<CvRect>& faceRect,std::string modelPath)
{
	int size;
	FILE* file;
	file = fopen(modelPath.c_str(), "rb");
	if (!file)
	{
		printf("# cannot read cascade from file %s",modelPath.c_str());
		return 1;
	}

	char input[1024], output[1024];

	fseek(file, 0L, SEEK_END);
	size = ftell(file);
	fseek(file, 0L, SEEK_SET);
	cascade = malloc(size);

	if (!cascade || size != fread(cascade, 1, size, file))
		return 1;

	fclose(file);

	//	set default parameters
	minsize = 25;
	maxsize = 30;

	angle = 0.0f;

	scalefactor = 1.3f;
	stridefactor = 0.1f;
	//modify by ruyiwei
	qthreshold = -99.0f;

	usepyr = 0;
	noclustering = 0;
	verbose = 0;

	//	//
	input[0] = 0;
	output[0] = 0;



	CvSize dst_cvsize;
	dst_cvsize.width = img->width;
	dst_cvsize.height = img->height;
	IplImage *dst = cvCreateImage(dst_cvsize, img->depth, img->nChannels);
	cvResize(img, dst, CV_INTER_LINEAR);
	
	if (!img)
	{
		printf("# cannot load image from image ");
		return 1;
	}
	try{
		clock_t t1, t2;
		t1 = clock();
		process_image(dst, faceRect, 1);
		t2 = clock();
		printf("# using time is '%d' ms\n",t2-t1);
	}
	catch (std::exception e){
		std::cout << e.what() <<std:: endl;
	}

	cvReleaseImage(&img);
	cvReleaseImage(&dst);
	return 0;
}
	




//int main(int argc, char* argv[])
//{
//	//
//	int arg;
//	char input[1024], output[1024];
//
//	{
//		int size;
//		FILE* file;
//
//		//
//		file = fopen("./resource/cascades/facefinder", "rb");
//
//		if (!file)
//		{
//			printf("# cannot read cascade from '%s'\n", argv[1]);
//			return 1;
//		}
//
//		//
//		fseek(file, 0L, SEEK_END);
//		size = ftell(file);
//		fseek(file, 0L, SEEK_SET);
//		//
//		cascade = malloc(size);
//
//		if (!cascade || size != fread(cascade, 1, size, file))
//			return 1;
//
//		//
//		fclose(file);
//	}
//
//	// set default parameters
//	minsize = 25;
//	maxsize = 30;
//
//	angle = 0.0f;
//
//	scalefactor = 1.3f;
//	stridefactor = 0.1f;
//
//	qthreshold = -99.0f;
//
//	usepyr = 0;
//	noclustering = 0;
//	verbose = 0;
//
//	//
//	input[0] = 0;
//	output[0] = 0;
//
//	// parse command line arguments
//	arg = 2;
//
//	if (verbose)
//	{
//		//
//		printf("# Copyright (c) 2017, RuMaxWell\n");
//		printf("# All rights reserved.\n\n");
//
//		printf("# cascade parameters:\n");
//		printf("#	tsr = %f\n", ((float*)cascade)[0]);//版本号
//		printf("#	tsc = %f\n", ((float*)cascade)[1]);//
//		printf("#	tdepth = %d\n", ((int*)cascade)[2]);
//		printf("#	ntrees = %d\n", ((int*)cascade)[3]);
//		printf("# detection parameters:\n");
//		printf("#	minsize = %d\n", minsize);
//		printf("#	maxsize = %d\n", maxsize);
//		printf("#	scalefactor = %f\n", scalefactor);
//		printf("#	stridefactor = %f\n", stridefactor);
//		printf("#	qthreshold = %f\n", qthreshold);
//		printf("#	usepyr = %d\n", usepyr);
//	}
//
//	{
//		IplImage* img;
//		cvNamedWindow("tmp",CV_WINDOW_KEEPRATIO);
//		
//		std::ifstream in(".//resource//Path_Images.txt");
//		std::string path;
//		while (in >> path){
//			std::cout << path << std::endl;
//			img = cvLoadImage(path.c_str(), CV_LOAD_IMAGE_COLOR);
//			CvSize dst_cvsize;
//			dst_cvsize.width = 256;
//			dst_cvsize.height = 256;
//			IplImage *dst = cvCreateImage(dst_cvsize, img->depth, img->nChannels);
//			cvResize(img, dst, CV_INTER_LINEAR);
//
//			if (!img)
//			{
//				printf("# cannot load image from '%s'\n", argv[3]);
//				return 1;
//			}
//			try{
//				clock_t t1, t2;
//				t1 = clock();
//				process_image(dst, 1);
//				t2 = clock();
//				printf("# using time is '%d' ms\n",t2-t1);
//			}
//			catch (std::exception e){
//				std::cout << e.what() <<std:: endl;
//			}
//			//
//			if (0 != output[0])
//				cvSaveImage(output, img, 0);
//			else if (!verbose)
//			{
//				cvShowImage("tmp", dst);
//				cvWaitKey(0);
//			}
//
//			//
//			cvReleaseImage(&img);
//			cvReleaseImage(&dst);
//		}
//	}
//
//	return 0;
//}
