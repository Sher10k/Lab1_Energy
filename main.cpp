#include <QCoreApplication>

#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

using namespace std;
using namespace cv;

// --
void erosion(const Mat &input_img, Mat &output_img, int apert)
{
    output_img = Mat::zeros(input_img.size(), CV_8U);
    for (int i = (apert / 2); i < input_img.cols - (apert / 2); i++) { 
        for (int j = (apert / 2); j < input_img.rows - (apert / 2); j++) {
            //uchar pix_value = input_img.at<uchar>(j, i);
            uchar min = 255;
            for (int ii = 0 - (apert / 2); ii <= (apert / 2); ii++) {
                for (int jj = 0 - (apert / 2); jj <= (apert / 2); jj++) {
                    uchar Y = input_img.at<uchar>(j + jj, i + ii);
                    if (Y < min) min = Y;
                }
            }
            output_img.at<uchar>(j, i) = min;
        }
    }
} 

// ++
void dilatac(const Mat &input_img, Mat &output_img, int apert) 
{
    output_img = Mat::zeros(input_img.size(), CV_8U);
    for (int i = (apert / 2); i < input_img.cols - (apert / 2); i++) {
        for (int j = (apert / 2); j < input_img.rows - (apert / 2); j++) {
            //uchar pix_value = input_img.at<uchar>(j, i);
            uchar max = 0;
            for (int ii = 0 - (apert / 2); ii <= (apert / 2); ii++) {
                for (int jj = 0 - (apert / 2); jj <= (apert / 2); jj++) {
                    uchar Y = input_img.at<uchar>(j + jj, i + ii);
                    if (Y > max) max = Y;
                }
            }
            output_img.at<uchar>(j, i) = max;
        } 
    } 
} 

void Roberts(const Mat &input_img, Mat &output_img)
{
    output_img = Mat::zeros(input_img.size(), CV_8U);
    float Rf1[3][3] = { {  1,  1,  1 },
                        {  1, -2,  1 },
                        { -1, -1, -1 } };
    
    float Rf2[3][3] = { { -1,  1,  1 },
                        { -1, -2,  1 },
                        { -1,  1,  1 } };
    
    float Rf3[3][3] = { { -1, -1,  1 },
                        { -1, -2,  1 },
                        {  1,  1,  1 } };
    for ( int i = 1; i < input_img.cols - 1; i++ ) 
    {
        for ( int j = 1; j < input_img.rows - 1; j++ ) 
        { 
            
            
            for ( int ii = -1; ii <= 1; ii++ )
            {
                for ( int jj = -1; jj <= 1; jj++ ) 
                {
                    
                }
            } 
        }
    }
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

        // Read image
    Mat img_in1 = imread( "image_1.png", 3 );
    Mat img_in2 = imread( "image_2.png", 3 );
    
        // Convert to gray
    Mat img_grey1, img_grey2;
    cvtColor( img_in1, img_grey1, COLOR_BGR2GRAY );
    cvtColor( img_in2, img_grey2, COLOR_BGR2GRAY );
    imwrite( "img_grey1.png", img_grey1 );
    imwrite( "img_grey2.png", img_grey2 );
    
        // Calculate difference between image
    Mat img_diff;
//    img_diff = img_grey1.clone();
//    for ( int i = 0; i < img_diff.total(); i++ )
//        img_diff.at< uchar >(i) = abs(img_grey2.at< uchar >(i) - img_grey1.at< uchar >(i));
    absdiff( img_grey1, img_grey2, img_diff );
    imwrite( "out_diff.png", img_diff );
    
        // Edge
    Mat img_edge1;
    Canny( img_grey1, img_edge1, 70, 150, 3, false );
    imwrite( "out_edge1.png", img_edge1 );
    
        // Binarization
    Mat img_edge1_bin, img_diff_bin;
    vector< int > params_write;
    params_write.push_back( IMWRITE_PNG_BILEVEL );
    threshold( img_edge1, img_edge1_bin, 17, 1, THRESH_BINARY  );   // THRESH_BINARY_INV THRESH_BINARY
    threshold( img_diff, img_diff_bin, 17, 1, THRESH_BINARY  );
    imwrite( "img_edge1_bin.png", img_edge1_bin, params_write );
    imwrite( "img_diff_bin.png", img_diff_bin, params_write );
    
        // Logical AND
    Mat img_and = img_diff.clone();
    for ( int i = 0; i < int(img_and.total()); i++ )
        img_and.at< uchar >(i) = img_diff_bin.at< uchar >(i) & img_edge1_bin.at< uchar >(i);
    imwrite( "img_and.png", img_and, params_write );
    
        // Morphological operations
    
    
    
    waitKey(0);
    return 0;   // a.exec();
}
