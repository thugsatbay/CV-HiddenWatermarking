//
// Watermark.cpp : Add watermark to an image, or inspect if a watermark is present.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// ctewani-gdhody-bansalro
//
//

//Link to the header file
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <SImage.h>
#include <SImageIO.h>
#include <fft.h>
#include <math.h>
#include <typeinfo>
#include <dirent.h>
#include <sys/types.h>


#define PI 3.14159265
#define WATERMARK_CONSTANT 0.25
#define WATERMARK_ALPHA 10
//RADIUS_OF_WATERMARK = (IMAGE_WIDTH/2) - WATERMARK_RADIUS_OFFSET
#define WATERMARK_RADIUS_OFFSET 50


using namespace std;

//Code Reference Started
// ---http://www.sanfoundry.com/c-program-integer-to-string-vice-versa/ ---
//Hacked By Gurleen Singh Dhody -gdhody- 12/Feb/2017
string toString(int num)
{
	int rem, len = 0, n;
    n = num;
    while (n != 0)
    {
        len++;
        n /= 10;
    }
    char *value = new char[len + 1];
    for (int i = 0; i < len; i++)
    {
        rem = num % 10;
        num = num / 10;
        value[len - i - 1] = rem + '0';
    }
    value[len] = '\0';
    //printf("Number is %s\n",value);
    string result(value);
    //printf("%s\n",result.c_str());
	return result;
}
//Code Reference Finished

// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const SDoublePlane &input, SDoublePlane &fft_real, SDoublePlane &fft_imag)
{
  fft_real = input;
  fft_imag = SDoublePlane(input.rows(), input.cols());

  FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const SDoublePlane &input_real, const SDoublePlane &input_imag, SDoublePlane &output_real)
{
  output_real = input_real;
  SDoublePlane output_imag = input_imag;

  FFT_2D(0, output_real, output_imag);
}


//Definition of function
SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag);


SDoublePlane check_image(const SDoublePlane &input, int N){
	srand(N);
	double *v,*c;
	SDoublePlane real,imagine;
	fft(input,real,imagine);
	// 4 is for quadrant division
	//---L---
	int l = N;
	v = new double[l];
	c = new double[l];
	for(int bitGen = 0;bitGen <= l; ++bitGen){
		v[bitGen] = (double) (rand() % 2);
	}
	//Constants Initialization
	int centerPoint = input.rows()/2 - 1;
	//---RADIUS---
	int radius = (int)(1 * (input.rows()/2)) - WATERMARK_RADIUS_OFFSET;
	//---ALPHA---
	double alpha = WATERMARK_ALPHA;
	//4 Quadrants
	int quadrantL = l / 4;
	//Angle, Real Value Holder, Absolute Constant
	double thetha = 0.0, Rvalue;
	int moveConstantOverV = 0;
	//Center of FFT
	int Xcenter = centerPoint, Ycenter = centerPoint;
	//To visualize circular watermark
	printf("Initialization For Watermark Done\n");
	//Finding Real values to find watermark
	for(int quad = 0; quad < 4; ++quad){
		for(int loop = 0; loop < quadrantL; ++loop){
			thetha = (PI* ((double)(quad*90.00) + ((loop+1)*(90.00/(double)quadrantL)) )) / 180.00;
			Rvalue = real[Xcenter + (int)ceil(radius * cos(thetha))][Ycenter + (int)ceil(radius * sin(thetha))];
			c[moveConstantOverV++] = Rvalue;
		}
	}

	//Calculating correlation between real values and vector v
	double c_mean, v_mean, c_std, v_std, cv_cov = 0.0;

	//Mean of c and v
	for(int loop = 0; loop < l; ++loop){
		c_mean += c[loop];
		v_mean += v[loop];
	}
	c_mean /= l;
	v_mean /= l;

	//Standard Deviation and covariance
	for(int loop = 0; loop < l; ++loop){
		c_std += (c[loop] - c_mean) * (c[loop] - c_mean);
		v_std += (v[loop] - v_mean) * (v[loop] - v_mean);
		cv_cov += (c[loop] - c_mean) * (v[loop] - v_mean);
	}

	//Pearson's Correlation r
	double pearsonCorrelation = cv_cov / sqrt(c_std * v_std);

	printf("Pearson's Correlation found to be with c and v vectors : %f\n",pearsonCorrelation);

	//Watermark present or not
	if (pearsonCorrelation >= WATERMARK_CONSTANT)
		printf("WaterMark Present in Image\n");
	else
		printf("Watermark Not Present in Image\n");

	return fft_magnitude(real,imagine);
}


SDoublePlane mark_image(const SDoublePlane &input, int N){
	//Generating the v vector through N and l
	srand(N);
	SDoublePlane finalC = SDoublePlane(input.rows(),input.cols());
	//---L---
	int l = N;
	double *v;
	v = new double[l];
	//Assigning the v vector with key N
	for(int bitGen = 0;bitGen <= l; ++bitGen){
		v[bitGen] = (double) (rand() % 2);
		//printf("%f",v[bitGen]);
	}
	//Constants Initialization
	SDoublePlane real,imagine;
	fft(input,real,imagine);
	int centerPoint = input.rows()/2 - 1;
	//---RADIUS---
	int radius = (int)(1 * (input.rows()/2)) - WATERMARK_RADIUS_OFFSET;
	//---ALPHA---
	double alpha = WATERMARK_ALPHA;
	//4 Quadrants
	int quadrantL = l / 4;
	//Angle, Real Value Holder, Absolute Constant
	double thetha = 0.0, Rvalue, negativeConstant = 1.0;
	int moveConstantOverV = 0;
	//Center of FFT
	int Xcenter = centerPoint, Ycenter = centerPoint;
	//To visualize circular watermark
	for(int i=0;i < input.rows();++i){
		for(int j=0;j < input.cols();++j){
			finalC[i][j] = 0;
		}
	}
	printf("Initialization For Watermark Done\n");
	//Updating Real Values to add watermark
	for(int quad = 0; quad < 4; ++quad){
		for(int loop = 0; loop < quadrantL; ++loop){
			//Calculating thetha to to find bin/point in our watermark circle to insert magnitude
			thetha = (PI* ((double)(quad*90.00) + ((loop+1)*(90.00/(double)quadrantL)) )) / 180.00;
			Rvalue = real[Xcenter + (int)ceil(radius * cos(thetha))][Ycenter + (int)ceil(radius * sin(thetha))];
			negativeConstant = 1.0;
			if (Rvalue < 0.0)
				negativeConstant = -1.0;
			//R(@) = R(@) + (alpha * |R(@)| * v[@])
			real[Xcenter + (int)ceil(radius * cos(thetha))][Ycenter + (int)ceil(radius * sin(thetha))] = Rvalue + (alpha * Rvalue * negativeConstant * v[moveConstantOverV]);
			//For Watermark_ILLUSTRATION
			finalC[Xcenter + (int)ceil(radius * cos(thetha))][Ycenter + (int)ceil(radius * sin(thetha))] = v[moveConstantOverV]*255;
			moveConstantOverV = moveConstantOverV + 1;
		}
	}

	//Printing Images FFT with watermark and its Illustration
	SDoublePlane watermarkFFT = fft_magnitude(real,imagine);
	SImageIO::write_png_file("Watermark_FFT.png",watermarkFFT,watermarkFFT,watermarkFFT);
	SDoublePlane outputReturnImage;
	ifft(real,imagine,outputReturnImage);
	SImageIO::write_png_file("Watermark_Circle_Illusration.png",finalC,finalC,finalC);
	printf("ADD-Watermark FFT and Circle Illustration Images Generated For Referral\n");
	return outputReturnImage;
}




SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag){

	//Declaring spectrogram that contains magnitude of FFT
	SDoublePlane spectrogram = SDoublePlane(fft_real.rows(),fft_real.cols());
  	//Declaring constants needed to run code
	int rowLoop = 0, columnLoop = 0;
	double minV = LONG_MAX,maxV = LONG_MIN, imageColorConstant = 255.00, value = 0.0;

	//Logic to calculate magnitude and find min and max value in our fft, min max used for normalization of image so we can display it correctly in output
	for (rowLoop = 0;rowLoop < fft_real.rows(); ++rowLoop){
		for(columnLoop = 0;columnLoop < fft_real.cols(); ++columnLoop){
			value = ( sqrt( (fft_real[rowLoop][columnLoop]*fft_real[rowLoop][columnLoop]) + (fft_imag[rowLoop][columnLoop]*fft_imag[rowLoop][columnLoop]) ) );

			//Log of 0 is -inf we don't want that
			if (value == 0)
				spectrogram[rowLoop][columnLoop] = LONG_MIN;
			else
				spectrogram[rowLoop][columnLoop] = log(value);

			//Find Min Max value except where value is LONG_MIN
			if (spectrogram[rowLoop][columnLoop] != LONG_MIN){
				if (spectrogram[rowLoop][columnLoop] < minV)
					minV = spectrogram[rowLoop][columnLoop];
				if (spectrogram[rowLoop][columnLoop] > maxV)
					maxV = spectrogram[rowLoop][columnLoop];
			}
		}
	}

	//updating the normalied values in spectrogram for displaying in png image
	for (rowLoop = 0;rowLoop < fft_real.rows(); ++rowLoop){
		for(columnLoop = 0;columnLoop < fft_real.cols(); ++columnLoop){
			//IF Long_Min means the magnitude was log of 0 so it is 0
			if (spectrogram[rowLoop][columnLoop] != LONG_MIN)
				spectrogram[rowLoop][columnLoop] = ((spectrogram[rowLoop][columnLoop] - minV) * (imageColorConstant/(maxV-minV)));
			else
				spectrogram[rowLoop][columnLoop] = 0;
		}
	}

	printf("FFT-Magnitude min max without normalization %f,%f; After normalization (0-255)\n",minV,maxV);
	printf("FFT-Magnitude Computed For Image FFT\n");
	return spectrogram;
}



SDoublePlane remove_interference(const SDoublePlane &input){
	SDoublePlane real,imagine,result,Mresult;
	//Do The FFT
	fft(input,real,imagine);
	int rowLoop = 0, columnLoop = 0;
	//Box Coordinates To Remove Message Put Into Image As Noise - HI
	int interferenceBox [] = {155,350};
	int rowWidth = 5, columnWidth = 8;
	//Box Coordinates Defined
		//Logic To Remove Box Coordinates with 0 intensity
		for (rowLoop = 0;rowLoop <= rowWidth; ++rowLoop){
			for (columnLoop = 0; columnLoop <= columnWidth; ++columnLoop){
				real[interferenceBox[0] + rowLoop][interferenceBox[0] + columnLoop] = 0.0;
				real[interferenceBox[1] + rowLoop + 1][interferenceBox[1] + columnLoop] = 0.0;
				imagine[interferenceBox[0] + rowLoop][interferenceBox[0] + columnLoop] = 0.0;
				imagine[interferenceBox[1] + rowLoop + 1][interferenceBox[1] + columnLoop] = 0.0;
			}
		}
	Mresult = fft_magnitude(real,imagine);
	printf("RMI-Interference Removed\n");
	SImageIO::write_png_file("Noise_Removed_fft.png",Mresult,Mresult,Mresult);
	printf("RMI-FFT Magnitude Image Generated After Noise Message - HI - Removal : %s\n","Noise_Removed_fft.png");
	//Parse The Image Back Through IFFT
	ifft(real,imagine,result);
	//Return The Grayscale Component To Be Written Back To PNG
	return result;
}

void quantitative_analysis(SDoublePlane input_image){
		for(int i = 1; i < 100; i++){
			// i is the value of N
			cout<<" ###### "<< i<<" ###### "<<endl;
			check_image(input_image, i);
		}
	}

int main(int argc, char **argv)
{
  try {

    if(argc < 4)
      {
	cout << "Insufficent number of arguments; correct usage:" << endl;
	cout << "    p2 problemID inputfile outputfile" << endl;
	return -1;
      }

    string part = argv[1];
    string inputFile = argv[2];
    string outputFile = argv[3];
    cout << "In: " << inputFile <<"  Out: " << outputFile << endl;
    SDoublePlane input_image;
    if(strcmp(argv[1], "1.0") != 0){
      printf("%s\n",inputFile.c_str());
      cout<<inputFile.c_str();
      input_image = SImageIO::read_png_file(inputFile.c_str());
    }

    if(strcmp(argv[1], "1.0") == 0){
			DIR   *d;
		   struct dirent *dir;
		   d = opendir("./Test_Images_2/");
		   string *file_names = new string[34];
		   int i = 0;
		   if (d)
		   {
		     while ((dir = readdir(d)) != NULL)
		     {
		       //if ((*dir->d_name).find(".png")!=-1){
		       printf("%s\n", dir->d_name);
		       file_names[i] = dir->d_name;
		       i++;
		     //}
		     }
		     cout<<i<<endl;
		     closedir(d);
		   }

      cout << "inside 1.0"<<endl;
      for(int i=0; i<34; i++){
        printf("%s\n", file_names[i].c_str());
        if(file_names[i].find(".png")!= -1){
              input_image = SImageIO::read_png_file(("Test_Images_2/" + file_names[i]).c_str());
              printf("hello\n");
              SDoublePlane Aresult = mark_image(input_image,90);
              SImageIO::write_png_file(("Test_Output/"+file_names[i]+toString(i)+".png").c_str(),Aresult,Aresult,Aresult);
              SDoublePlane Aresult_check = check_image(input_image,atoi(argv[5]));
              SImageIO::write_png_file(("Test_Output/"+file_names[i]+toString(i)+"result.png").c_str(),Aresult_check,Aresult_check,Aresult_check);
      }
        //printf ("Watermark Add 1.3 Module Finished, with output Image : %s\n",argv[3]);
      }
    }
    else if(strcmp(argv[1],"1.1")==0)
      {
				SDoublePlane real,imagine,result;
				//Do the fft
				fft(input_image,real,imagine);
				//Find the magnitude
				result = fft_magnitude(real,imagine);
				SImageIO::write_png_file(argv[3],result,result,result);
				printf ("Magnitude 1.1 Module Finished, with output Image : %s\n",argv[3]);
      }
    else if(strcmp(argv[1],"1.2")==0)
      {
				SDoublePlane Nresult;
				//Remove the noise
				Nresult = remove_interference(input_image);
				//Save The Image
				SImageIO::write_png_file(argv[3],Nresult,Nresult,Nresult);
				printf ("Noise Removal 1.2 Module Finished, with output Image : %s\n",argv[3]);
      }
    else if(part == "1.3")
      {
	if(argc < 6)
	  {
	    cout << "Need 6 parameters for watermark part:" << endl;
	    cout << "    p2 1.3 inputfile outputfile operation N" << endl;
	    return -1;
	  }
	if(strcmp(argv[4],"add")==0)
	  {
		SDoublePlane Aresult = mark_image(input_image,atoi(argv[5]));
		SImageIO::write_png_file(argv[3],Aresult,Aresult,Aresult);
		printf ("Watermark Add 1.3 Module Finished, with output Image : %s\n",argv[3]);
	  }
	else if(strcmp(argv[4],"check")==0)
	  {
		SDoublePlane Aresult = check_image(input_image,atoi(argv[5]));
		SImageIO::write_png_file(argv[3],Aresult,Aresult,Aresult);
		printf ("Watermark Check 1.4 Module Finished, with output Image : %s\n",argv[3]);
    	  }
	else if(strcmp(argv[4],"quantCheck")==0)
	  {
		SDoublePlane result = mark_image(input_image, atoi(argv[5]));
		SImageIO::write_png_file(argv[3],result,result,result);

		quantitative_analysis(SImageIO::read_png_file(argv[3]));
  }
	else
	  throw string("Bad operation!");
      }
    else
      throw string("Bad part!");

  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}
