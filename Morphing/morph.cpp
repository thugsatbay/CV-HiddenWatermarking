#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <SImage.h>
#include <SImageIO.h>
#include <fft.h>

using namespace std;

int main(int argc, char** argv){

	//Gaussian Filter 5*5 Sigma = 10
	SDoublePlane gaussian = SDoublePlane(5,5);
	gaussian[0][0] = 0.039206;
	gaussian[0][1] = 0.039798;
	gaussian[0][2] = 0.039997;
	gaussian[0][3] = 0.039798;
	gaussian[0][4] = 0.039206;
	gaussian[1][0] = 0.039798;
	gaussian[1][1] = 0.040399;
	gaussian[1][2] = 0.040601;
	gaussian[1][3] = 0.040399;
	gaussian[1][4] = 0.039798;
	gaussian[2][0] = 0.039997;
	gaussian[2][1] = 0.040601;
	gaussian[2][2] = 0.040804;
	gaussian[2][3] = 0.040601;
	gaussian[2][4] = 0.039997;
	gaussian[3][0] = 0.039798;
	gaussian[3][1] = 0.040399;
	gaussian[3][2] = 0.040601;
	gaussian[3][3] = 0.040399;
	gaussian[3][4] = 0.039798;
	gaussian[4][0] = 0.039206;
	gaussian[4][1] = 0.039798;
	gaussian[4][2] = 0.039997;
	gaussian[4][3] = 0.039798;
	gaussian[4][4] = 0.039206;
	
	//Two Input Files Names To Be Morphed
	string inputFile_1 = argv[1];
	string inputFile_2 = argv[2];

	//Double Dimension Array To Process Images
	SDoublePlane input_1 = SImageIO::read_png_file(inputFile_1.c_str());
	SDoublePlane input_2 = SImageIO::read_png_file(inputFile_2.c_str());
	SDoublePlane lowPass_1 = SDoublePlane(input_1.rows(), input_1.cols());
	SDoublePlane lowPass_2 = SDoublePlane(input_1.rows(), input_1.cols());
	SDoublePlane highPass = SDoublePlane(input_1.rows(), input_1.cols());
	SDoublePlane morphPass = SDoublePlane(input_1.rows(), input_1.cols());

	//Convolution Gaussian Procedure - Low Pass Filter
	for(int i=2 ; i<input_1.rows()-2; i++){
		for(int j=2; j<input_1.cols()-2; j++){
			for(int k=0; k<gaussian.rows(); k++){
				for(int l=0; l<gaussian.cols(); l++){
					lowPass_1[i][j] = lowPass_1[i][j] + input_1[-2+k+i][-2+l+j] * (double)gaussian[gaussian.rows()-k-1][gaussian.cols()-l-1];
				}
			}	
		}
	}
	//Low Pass 1st Image
	SImageIO::write_png_file("lowPass_1.png", lowPass_1, lowPass_1, lowPass_1);

	//Convolution Gaussian Procedure - Low Pass Filter
	for(int i=2 ; i<input_2.rows()-2; i++){
		for(int j=2; j<input_2.cols()-2; j++){
			for(int k=0; k<gaussian.rows(); k++){
				for(int l=0; l<gaussian.cols(); l++){
					lowPass_2[i][j] += input_2[-2+i+k][-2+j+l] * (double)gaussian[gaussian.rows()-k-1][gaussian.cols()-l-1];
				}
			}	
		}
	}
	//Low Pass 2nd Image
	SImageIO::write_png_file("lowPass_2.png", lowPass_2, lowPass_2, lowPass_2);

	//High Pass = OriginalImage - LowImage
	for(int i=1; i<input_2.rows()-1; i++){
		for(int j=1; j<input_2.cols()-1; j++){
			highPass[i][j] = input_2[i][j] - lowPass_2[i][j];
		}
	}
	//High Pass Image
	SImageIO::write_png_file("highPass.png", highPass, highPass, highPass);

	//Morph Image = LowPass + HighPass
	for(int i=1; i<input_2.rows()-1; i++){
		for(int j=1; j<input_2.cols()-1; j++){
			morphPass[i][j] = (lowPass_1[i][j] + highPass[i][j]) / 2;
		}
	}
	
	//Output Image
	SImageIO::write_png_file("output.png", morphPass, morphPass, morphPass);
	return 0;
}
