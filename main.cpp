#include "faceDetectLandmark.h"

int main(int argc, char **argv) {


	fstream file;
	file.open("./list.txt");
	string modelPath = "./resource/cascades/facefinder";
 
	while (!file.eof())
	{
		cout<<"begin to detect"<<endl;
    string imageDir;
		file >> imageDir; 
    cout<<imageDir<<endl;
		Mat image = imread(imageDir, -1);
		if (image.empty()) continue;
		vector<Mat > landmarkPoint;
		vector<CvRect > rects;
		resize(image, image, Size(160, 240));
		run(image, rects, landmarkPoint,modelPath);
	}
	
	return 0;
}


//int main(int argc, char **argv) {
//    if (argc != 2) {
//        LOG("We need an argument");
//        return 0;
//    }
//    if (strcmp(argv[1], "train") == 0) {
//        return train(0);
//    }
//    if (strcmp(argv[1], "resume") == 0) {
//        int start_from;
//        printf("Which stage you want to resume from: ");
//        scanf("%d", &start_from);
//        return train(start_from);
//    }
//    if (strcmp(argv[1], "test") == 0) {
//        return test();
//    }
//    if (strcmp(argv[1], "prepare") == 0) {
//        return prepare();
//    }
//    if (strcmp(argv[1], "run") == 0) {
//		fstream file;
//		file.open("list.txt");
//		while (!file.eof())
//		{
//			string imageDir;
//			file >> imageDir;
//			Mat image = imread(imageDir, -1);
//			vector<int> landmarkPoint;
//			Rect rects;
//			run(image, rects, landmarkPoint);
//		}
//    }
//    else {
//        LOG("Wrong Arguments.");
//    }
//    return 0;
//}
