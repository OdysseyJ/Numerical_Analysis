//
//  main.cpp
//  opencv
//
//  Created by MacBookPro on 19/11/2019.
//  Copyright © 2019 MacBookPro. All rights reserved.
//
#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

int main(int argc, const char * argv[]) {
    
    // 크기 입력받기.
    double size;
    printf("Input the size: ");
    scanf("%lf", &size);
    
    // 원본 tiger.jpg 이미지.
    Mat image;
    image = imread("/Users/macbook/Desktop/opencv/opencv/tiger.jpg", IMREAD_COLOR);
    
    // 이미지 못찾은 경우 예외처리
    if (image.empty())
    {
        cout << "Could not open or find the image" << endl;
        return -1;
    }
    
    // 출력하고자 하는 이미지 파일
    Mat scaledImage(image.rows * size, image.cols * size, CV_8UC3, Scalar(0));
    
    // image scaling
    for(int r = 0; r < scaledImage.rows; r++)
    {
        for(int c = 0; c < scaledImage.cols; c++)
        {
            int pc = (int)(c / size);
            int pr = (int)(r / size);

            double fc1 = (double)c / (double)size - (double)pc;
            double fc2 = 1 - fc1;
            double fr1 = (double)r / (double)size - (double)pr;
            double fr2 = 1 - fr1;

            double w1 = fc2 * fr2;
            double w2 = fc1 * fr2;
            double w3 = fc2 * fr1;
            double w4 = fc1 * fr1;

            Vec3b P1 = image.at<Vec3b>(pr, pc);
            Vec3b P2 = image.at<Vec3b>(pr, pc < image.cols - 1? pc + 1 : pc);
            Vec3b P3 = image.at<Vec3b>(pr < image.rows - 1? pr + 1 : pr, pc);
            Vec3b P4 = image.at<Vec3b>(pr < image.rows - 1? pr + 1 : pr, pc < image.cols - 1? pc + 1 : pc);
            scaledImage.at<Vec3b>(r, c) = w1 * P1 + w2 * P2 + w3 * P3 + w4 * P4;
        }
    }
    
    // 이미지 해당 경로에 만들기.
    imwrite("/Users/macbook/Desktop/opencv/opencv/scaledTiger.jpg", scaledImage);

    return 0;
}
