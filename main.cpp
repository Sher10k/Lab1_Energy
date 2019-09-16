#include <QCoreApplication>

#include <vector>
#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

using namespace std;
using namespace cv;

// --
void erosion(const Mat &input_img, Mat &output_img, int apert)
{
    output_img = Mat::zeros( input_img.size(), input_img.type() );
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
    output_img = Mat::zeros( input_img.size(), input_img.type() );
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
    output_img = Mat::zeros( input_img.size(), input_img.type() );
    Mat Rf1 = ( Mat_< char >(3,3) << 1,  1, 1, 
                                      1, -2, 1, 
                                     -1, -1,-1 );
    Mat Rf2 = ( Mat_< char >(3,3) << -1,  1, 1, 
                                      -1, -2, 1, 
                                      -1,  1, 1 );
    Mat Rf3 = ( Mat_< char >(3,3) << -1, -1, 1, 
                                      -1, -2, 1, 
                                       1,  1, 1 );
    
    Mat img_Rf1, img_Rf2, img_Rf3;
    
    filter2D( input_img, img_Rf1, -1, Rf1 );
    filter2D( input_img, img_Rf2, -1, Rf2 );
    filter2D( input_img, img_Rf3, -1, Rf3 );
    
    addWeighted( img_Rf1, 1.0, img_Rf2, 1.0, 0, output_img );
    addWeighted( output_img, 1.0, img_Rf3, 1.0, 0, output_img );
}

void Sobel(const Mat &input_img, Mat &output_img)
{
    output_img = Mat::zeros( input_img.size(), input_img.type() );   
    int Fk[3][3] = { { 1,1,1 },
                     { 1,1,1 },    
                     { 1,1,1 } }; // маска фильтра 3*3 
    
    for (int i = 1; i < input_img.cols - 1; i++)
        for (int j = 1; j < input_img.rows - 1; j++) {
            // далее производим свертку
            for (int ii = -1; ii <= 1; ii++)
                for (int jj = -1; jj <= 1; jj++) {
                    //uchar blurred = input_img.at<uchar>(j + jj, i + ii);
                    Fk[ii + 1][jj + 1] = input_img.at< uchar >(j + jj, i + ii);
                }
            int blurred = abs((Fk[0][0] + 2 * Fk[0][1] + Fk[0][2]) - (Fk[2][0] + 2 * Fk[2][1] + Fk[2][2])) + 
                            abs((Fk[0][0] + 2 * Fk[1][0] + Fk[2][0]) - (Fk[0][2] + 2 * Fk[1][2] + Fk[2][2]));
            if (blurred < 0) blurred = 0;
            else if (blurred > 255) blurred = 255;
            //else blurred = blurred;
            output_img.at< uchar >(j, i) = uchar(blurred);
        }
} 

void Hist( Mat img, Mat &hh, Mat &hv)
{
    hh = Mat::ones( img.size(), img.type() );
    hv = Mat::ones( img.size(), img.type() );
    hh *= 255;
    hv *= 255;

    bitwise_not( img, img );
    
    int intensity = 0;
    for ( int i = 0; i < img.cols; i++ )
    {
        intensity = countNonZero( img.col(i) );
        for ( int k = intensity; k != 0; k-- )
            hh.at< uchar >(hh.rows - k, i) = 0;
    }
    for ( int j = 0; j < img.rows; j++ )
    {
        intensity = countNonZero( img.row(j) );
        for ( int k = intensity; k != 0; k-- )
            hv.at< uchar >(j, hv.cols - k) = 0;
    }
}

Mat draw_arctangle_new( const Mat img_horizon, const Mat img_vertical, Mat img )
{
	int level1 = 1;
	int level2 = 1;
    vector< int > X;
    vector< int > Y;
    
    bitwise_not( img_horizon, img_horizon );
    bitwise_not( img_vertical, img_vertical );
    bool truelevel = false;
    for ( int i = 0; i < img.cols; i++ )
    {
        int intensity = countNonZero( img_horizon.col(i) );
        if ( (intensity >= level1) && (!truelevel) ) 
        {
            X.push_back( i );
            truelevel = true;
        }
        else if ( (intensity < level1) && (truelevel) )
        {
            X.push_back( i-1 );
            truelevel = false;
        }
        if ( (intensity >= level1) && (i == img.cols-1) )
            X.push_back( i );
    }
    truelevel = false;
    for ( int j = 0; j < img.rows; j++ )
    {
        int intensity = countNonZero( img_vertical.row(j) );
        if ( (intensity >= level2) && (!truelevel) ) 
        {
            Y.push_back( j );
            truelevel = true;
        }
        else if ( (intensity < level2) && (truelevel) )
        {
            Y.push_back( j-1 );
            truelevel = false;
        }
        if ( (intensity >= level2) && (j == img.rows-1) )
            Y.push_back( j );
    }
    rectangle(img, Point(X.at(0), Y.at(0)), Point(X.at(1), Y.at(1)), Scalar(0, 0, 255), 2, 8, 0);
    rectangle(img, Point(X.at(2), Y.at(0)), Point(X.at(3), Y.at(1)), Scalar(0, 0, 255), 2, 8, 0);
	return img;
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

        // Read image
    Mat img_in1 = imread( "image_1.png", 3 );
    Mat img_in2 = imread( "image_2.png", 3 );
    
        // Crop
    img_in1 = img_in1( Rect( 200, 100, img_in1.cols - 400, img_in1.rows - 200 ) );
    img_in2 = img_in2( Rect( 200, 100, img_in2.cols - 400, img_in2.rows - 200 ) );
    
        // Convert to gray
    Mat img_grey1, img_grey2;
    cvtColor( img_in1, img_grey1, COLOR_BGR2GRAY );
    cvtColor( img_in2, img_grey2, COLOR_BGR2GRAY );
    imwrite( "img_grey1.png", img_grey1 );
    imwrite( "img_grey2.png", img_grey2 );
    
        // Gauss
    int sigma = 3;
    int ksize = 5;  //( sigma*5 ) | 1;
    GaussianBlur( img_grey1, img_grey1, Size( ksize, ksize ), sigma, sigma, cv::BORDER_DEFAULT );
    GaussianBlur( img_grey2, img_grey2, Size( ksize, ksize ), sigma, sigma, cv::BORDER_DEFAULT );
    imwrite( "img_gauss1.png", img_grey1 );
    imwrite( "img_gauss2.png", img_grey2 );
    
        // Calculate difference between image
    Mat img_diff;
//    img_diff = img_grey1.clone();
//    for ( int i = 0; i < img_diff.total(); i++ )
//        img_diff.at< uchar >(i) = abs(img_grey2.at< uchar >(i) - img_grey1.at< uchar >(i));
    absdiff( img_grey1, img_grey2, img_diff );
    imwrite( "img_diff.png", img_diff );
    
        // Edge
    Mat img_edge1;
//    Canny( img_grey1, img_edge1, 70, 150, 3, false );
//    Sobel( img_grey1, img_edge1 );
    Roberts( img_grey1, img_edge1 );
    imwrite( "img_edge1.png", img_edge1 );
    
        // Binarization
    Mat img_edge1_bin, img_diff_bin;
    vector< int > params_write;
    params_write.push_back( IMWRITE_PNG_BILEVEL );
    
    threshold( img_diff, img_diff_bin, 10, 255, THRESH_BINARY  );
    Mat idbe;
    erosion( img_diff_bin, idbe, 3 );
    idbe.copyTo( img_diff_bin );
    imwrite( "img_diff_bin.png", img_diff_bin );
    
    threshold( img_edge1, img_edge1_bin, 10, 255, THRESH_BINARY  );   // THRESH_BINARY_INV THRESH_BINARY
    imwrite( "img_edge1_bin.png", img_edge1_bin );
    
        // Logical AND
    Mat img_and = img_diff.clone();
//    for ( int i = 0; i < int(img_and.total()); i++ ) 
//        img_and.at< uchar >(i) = img_diff_bin.at< uchar >(i) & img_edge1_bin.at< uchar >(i);
    bitwise_and( img_diff_bin, img_edge1_bin, img_and );
    bitwise_not( img_and, img_and );    // invertion image
    imwrite( "img_and.png", img_and );  // , params_write
    
        // Morphological operations
    Mat img_temp = Mat::zeros( img_and.size(), img_and.type() );
    Mat img_morf = Mat::zeros( img_and.size(), img_and.type() );
    Mat element = getStructuringElement( MORPH_RECT, Size(10, 10), Point(-1, -1) );
    morphologyEx( img_and, img_temp, MORPH_OPEN, element );
    morphologyEx( img_temp, img_morf, MORPH_CLOSE, element );
    imwrite( "img_morfolog_1.png", img_morf );
//    morphologyEx( img_morf, img_temp, MORPH_OPEN, element );
//    morphologyEx( img_temp, img_morf, MORPH_CLOSE, element );
//    imwrite( "img_morfolog_2.png", img_morf );
    
        // Histogram
    Mat hist_vertical, hist_horizon;
    Hist( img_morf, hist_horizon, hist_vertical );
    imwrite( "hist_horizon.png", hist_horizon );
    imwrite( "hist_vertical.png", hist_vertical );
    
        // Strobe
    Mat img_finish;
    img_finish = draw_arctangle_new( hist_horizon, hist_vertical, img_in1 );
    imwrite( "finish.png", img_finish );
    
    waitKey(0);
    return 0;   // a.exec();
}


