#define _USE_MATH_DEFINES
#include <opencv2/opencv.hpp>
#include <iostream>
#include <ctime>

using namespace cv;
using namespace std;

const float pi = float(M_PI);

Mat DftFilter(int rows, int cols, float sigma)
{
	int m = getOptimalDFTSize(rows);
	int n = getOptimalDFTSize(cols);

	Mat h(m, n, CV_32F);
	for (int i = 0; i < m; i++) {
		auto hrow = h.ptr<float>(i);
		for (int j = 0; j < n; j++) {
			float ux = j < n / 2 ? 0.f : float(n - 1);
			float uy = i < m / 2 ? 0.f : float(m - 1);
			float val = exp(-2 * pi*pi*sigma*sigma*((j-ux)*(j-ux) + (i-uy)*(i-uy)) / (m*n));
			hrow[j] = val;
		}
	}
	return h;
}


Mat DftBlur(Mat src, Mat filter)
{
	int m = getOptimalDFTSize(src.rows);
	int n = getOptimalDFTSize(src.cols);

	assert(filter.rows == m);
	assert(filter.cols == n);

	int k = src.channels();
	vector<Mat> rgb(k);
	split(src, rgb);
	for (int i = 0; i < k; i++) {
		Mat I = rgb[i];
		copyMakeBorder(I, I, 0, m - src.rows, 0, n - src.cols, BORDER_CONSTANT, Scalar::all(0));
		Mat planes[] = { I, Mat::zeros(I.size(), CV_32F) };
		Mat F;
		merge(planes, 2, F);
		dft(F, F);
		split(F, planes);
		Mat disp;
		log(planes[0], disp);
		disp = disp / 10.0;

		planes[0] = planes[0].mul(filter);
		planes[1] = planes[1].mul(filter);
		merge(planes, 2, F);
		cv::dft(F, I, DFT_INVERSE + DFT_REAL_OUTPUT + DFT_SCALE);
		rgb[i] = I;
	}
	Mat dst;
	merge(rgb, dst);
	Rect roi(0, 0, src.cols, src.rows);
	return dst(roi);
}


Mat DctFilter(int rows, int cols, float sigma)
{
	int m = 2 * getOptimalDFTSize((rows + 1) / 2);
	int n = 2 * getOptimalDFTSize((cols + 1) / 2);

	Mat h(m, n, CV_32F);
	for (int i = 0; i < m; i++) {
		auto hrow = h.ptr<float>(i);
		for (int j = 0; j < n; j++) {
			float val = exp(-0.5f * pi*pi*sigma*sigma*(i*i + j*j) / (m*n));
			hrow[j] = val;
		}
	}
	return h;
}

Mat DctBlur(Mat src, Mat filter)
{
	int m = 2 * getOptimalDFTSize((src.rows + 1) / 2);
	int n = 2 * getOptimalDFTSize((src.cols + 1) / 2);

	assert(filter.rows == m);
	assert(filter.cols == n);

	int k = src.channels();
	vector<Mat> rgb(k);
	split(src, rgb);
	for (int i = 0; i < k; i++) {
		Mat C = rgb[i];
		copyMakeBorder(C, C, 0, m - src.rows, 0, n - src.cols, BORDER_REPLICATE);
		dct(C, C);
		C = C.mul(filter);
		idct(C, C);
		rgb[i] = C;
	}
	Mat dst;
	merge(rgb, dst);
	Rect roi(0, 0, src.cols, src.rows);
	return dst(roi);
}

int main()
{
	Mat Iin = imread("input.jpg");
	Iin.convertTo(Iin, CV_32FC3, 1 / 255.0);
	cvtColor(Iin, Iin, CV_BGR2RGB);

	int K = 16;
	for (int i = 0; i < K; i++) {
		float sigma = 10;// pow(2.f, i - 1);
		float ratio = pow(2.f, 0.2f*(i - 8));
		int rows = int(Iin.rows*ratio);
		int cols = int(Iin.cols*ratio);
		Mat src;
		resize(Iin, src, Size(cols, rows));
		//cout << sigma << "\t";
		cout << rows << '\t' << cols << '\t' << rows*cols << '\t';

		clock_t s, e;
		s = clock();
		Mat spatial_dst;
		GaussianBlur(src, spatial_dst, Size(), sigma);
		e = clock();
		cout << float(e - s) * 1000 / CLOCKS_PER_SEC << '\t';

		//imshow("spatial", spatial_dst);

		Mat dft_filter = DftFilter(src.rows, src.cols, sigma);
		s = clock();
		Mat dft_dst = DftBlur(src, dft_filter);
		e = clock();
		cout << float(e - s) * 1000 / CLOCKS_PER_SEC << '\t';
		//imshow("dft", dft_dst);

		Mat dct_filter = DctFilter(src.rows, src.cols, sigma);
		s = clock();
		Mat dct_dst = DctBlur(src, dct_filter);
		e = clock();
		cout <<  float(e - s) * 1000 / CLOCKS_PER_SEC << '\t' ;
		//imshow("dct", dct_dst);
		//waitKey();
		cout << ';' << endl;
	}
	system("pause");

	return 0;
}