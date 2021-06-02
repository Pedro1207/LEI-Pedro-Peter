//
// Lighthouse3D.com VS*L OpenGL Sample
//
// Loading and displaying a Textured Model
//
// Uses:
//  Assimp 3.0 library for model loading
//		http://assimp.sourceforge.net/
//  Devil for image loading
//		http://openil.sourceforge.net/
//  GLEW for OpenGL post 1.1 functions
//		http://glew.sourceforge.net/
//	TinyXML for font definition parsing
//		http://sourceforge.net/projects/tinyxml/
//
// This demo was built for learning purposes only.
// Some code could be severely optimised, but I tried to
// keep as simple and clear as possible.
//
// The code comes with no warranties, use it at your own risk.
// You may use it, or parts of it, wherever you want.
//
// If you do use it I would love to hear about it. 
// Just post a comment at Lighthouse3D.com

// Have Fun :-)

#include <math.h>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define M_PIf float(M_PI)


// include GLEW to access OpenGL 3.3 functions
#include <GL/glew.h>

// GLUT is the toolkit to interface with the OS
#include <GL/freeglut.h>

// Use Very Simple Libs
#include <vsl/vslibs.h>

#include <fftw/fftw3.h>

#include "config.h"

VSMathLib *vsml;
VSShaderLib program, programFonts, programFFT_H, programFFT_V, programHTK;

#if (__VSL_FONT_LOADING__ == 1) && (__VSL_TEXTURE_LOADING__ == 1)
VSFontLib vsfl;
#endif

VSModelLib myModel;

VSAxis axis;
VSGrid gridY;

unsigned int aSentence, profileSentence;

// Query to track the number of primitives
// issued by the tesselation shaders
GLuint counterQ;
unsigned int primitiveCounter = 0;

// Camera Position
float camX = 0, camY = 0, camZ = 5;

// Mouse Tracking Variables
int startX, startY, tracking = 0;

// Camera Spherical Coordinates
float alpha = 0.0f, beta = 0.0f;
float r = 50.0f;

//// Frame counting and FPS computation
long myTime,timebase = 0,frame = 0;
char s[256];

float lightDir[4] = { 1.0f, 1.0f, 1.0f, 0.0f };




struct point2d {

	float x, y;

	point2d() { x = 0.0f; y = 0.0f; }
	point2d(float x, float y) { this->x = x; this->y = y; }

	float length() {
		return sqrt(x*x + y*y);
	}

	point2d normalize() {

		float n = sqrt(x*x + y*y);
		if (n != 0)
			return point2d(x / n, y / n);
		else
			return point2d(0, 0);
	}

	float dot(point2d &a) {
		return a.x*x + a.y*y;
	}

	point2d operator-() const {
		return point2d(-x, -y);
	}
};


struct complex {

	float r, i;

	complex() { r = 0.0f; i = 0.0f; }
	complex(point2d &p) { r = p.x ; i = p.y; }
	complex(float r, float i) { this->r = r; this->i = i; }

	float mag() {
		return (sqrt(r*r + i*i));
	}

	float phase() {
		return atan(i / r);
	}

	complex operator +(const complex& a) const {
		return complex(r + a.r, i + a.i);
	}

	complex operator -(const complex& a) const {
		return complex(r - a.r, i - a.i);
	}

	complex operator *(const complex& a) const {
		return complex(a.r * r - a.i * i, a.r * i + a.i * r);
	}

	// c * k
	complex operator *(const float a) const {
		return complex(a * r, a * i);
	}
};


// allows for k * c;
complex operator *(float a, const complex& c) {
	return complex(a * c.r, a * c.i);
}

// returns e^ki
complex euler(float k) {
	return complex(cosf(k), sinf(k));
}

complex conjugate(const complex& c) {
	return complex(c.r, -c.i);
}

struct point3d {
	float x, y, z;
};

int L = 1024;
#define FW 256// Fourier dimension FW x FW

float max, min;

// Chen, Li-ning and Jin, Yi - cheng and Yin, Yong and Ren, Hong-xiang
// On the Wave Spectrum Selection in Ocean Wave Scene Simulation of the Maritime Simulator
// AsiaSim 2013
#define A 3.48e-3f


struct point2d windDir = point2d(1.0f, 0.0f);
float windSpeed = 25;

complex h0k[FW * FW];
complex h0kconj[FW * FW];

complex htk[FW * FW];
complex slopeX[FW * FW], slopeZ[FW * FW];
complex dX[FW * FW], dZ[FW * FW];

float pos[FW * FW * 4];
float normal[FW * FW * 3];
float tc[FW * FW * 2];



complex *outH, *outSX, *outSZ, *outDX, *outDZ ;

VSSolidGrid sg;
float timer = 0000.0f; // time

int useFFTW = 3;
int htkCPU = 1;



// -------------------------------------------------------------------------

void displayArray(complex *c, int nx, int ny) {

	for (int i = 0; i < nx; ++i) {

		for (int j = 0; j < ny; ++j) {
			printf("%f - %f    ", c[i*ny + j].r, c[i*ny + j].i);
		}
		printf("\n");
	}
	printf("\n");
}


/*-------------------------------------------------------------------------
http://paulbourke.net/miscellaneous/dft/
Perform a 2D FFT inplace given a complex 2D array
The direction dir, 1 for forward, -1 for reverse
The size of the array (nx,ny)
Return false if there are memory problems or
the dimensions are not powers of 2
*/
int FFT(int dir, int m, double *x, double *y);


int FFT2D(complex *c, int nx, int ny, int dir)
{
	int i, j;
	int m = log2(FW), twopm;
	double *real, *imag;

//	displayArray(c, nx, ny);


	/* Transform the columns */
	real = (double *)malloc(ny * sizeof(double));
	imag = (double *)malloc(ny * sizeof(double));
	for (i = 0; i<nx; i++) {
		for (j = 0; j<ny; j++) {
			real[j] = c[i*FW + j].r;
			imag[j] = c[i*FW + j].i;
		}
		FFT(dir, m, real, imag);
		for (j = 0; j<ny; j++) {
			c[i*FW + j].r = real[j];
			c[i*FW + j].i = imag[j];
		}
	}
	free(real);
	free(imag);
	//displayArray(c, nx, ny);

	/* Transform the rows */
	real = (double *)malloc(nx * sizeof(double));
	imag = (double *)malloc(nx * sizeof(double));
	for (j = 0; j<ny; j++) {
		for (i = 0; i<nx; i++) {
			real[i] = c[i*FW + j].r;
			imag[i] = c[i*FW + j].i;
		}
		FFT(dir, m, real, imag);
		for (i = 0; i<nx; i++) {
			c[i*FW + j].r = real[i];
			c[i*FW + j].i = imag[i];
		}
	}
	free(real);
	free(imag);

	//displayArray(c, nx, ny);

	return(TRUE);
}

/*-------------------------------------------------------------------------
This computes an in-place complex-to-complex FFT
x and y are the real and imaginary arrays of 2^m points.
dir =  1 gives forward transform
dir = -1 gives reverse transform

Formula: forward
N-1
---
1   \          - j k 2 pi n / N
X(n) = ---   >   x(k) e                    = forward transform
N   /                                n=0..N-1
---
k=0

Formula: reverse
N-1
---
\          j k 2 pi n / N
X(n) =       >   x(k) e                    = forward transform
/                                n=0..N-1
---
k=0
*/
int FFT(int dir, int m, double *x, double *y)
{
	long nn, i, i1, j, k, i2, l, l1, l2;
	double c1, c2, tx, ty, t1, t2, u1, u2, z;

	/* Calculate the number of points */
	nn = 1;
	for (i = 0; i<m; i++)
		nn *= 2;

	/* Do the bit reversal */
	i2 = nn >> 1;
	j = 0;
	for (i = 0; i<nn - 1; i++) {
		if (i < j) {
			tx = x[i];
			ty = y[i];
			x[i] = x[j];
			y[i] = y[j];
			x[j] = tx;
			y[j] = ty;
			//printf("%d %d\n", i, j);
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	/* Compute the FFT */
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	for (l = 0; l<m; l++) {
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;
		for (j = 0; j<l1; j++) {
			for (i = j; i<nn; i += l2) {
				i1 = i + l1;
				t1 = u1 * x[i1] - u2 * y[i1];
				t2 = u1 * y[i1] + u2 * x[i1];
				//printf("%d %d - %f %f - %f %f - %f %f\n", i, i1, x[i], y[i], x[i1], y[i1], u1, u2);
				x[i1] = x[i] - t1;
				y[i1] = y[i] - t2;
				x[i] += t1;
				y[i] += t2;
				//printf("%d %d - %f %f - %f %f\n", i, i1, x[i], y[i], x[i1], y[i1]);
			}
			z = u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		if (dir == 1)
			c2 = -c2;
		c1 = sqrt((1.0 + c1) / 2.0);

		//for (int i = 0; i < 8; ++i) {
		//	printf("#%f - %f    ", x[i], y[i]);
		//}
		//printf("\n\n");
	}

	/* Scaling for forward transform */
	if (dir == 1) {
		for (i = 0; i<nn; i++) {
			x[i] /= (double)nn;
			y[i] /= (double)nn;
		}
	}

	return(TRUE);
}






// ------------------------------------------------------------
//
// My FFT2D
//

// bit reverse from https://alikhuram.wordpress.com/2013/05/15/computing-the-fast-fourier-transform-algorithm-in-c/

std::vector<int> br;
complex preW[FW][FW+1];


unsigned int bitReverse(unsigned int x, int log2n)

{
	int n = 0;
	int mask = 0x1;
	for (int i = 0; i < log2n; i++)

	{
		n <<= 1;
		n |= (x & 1);
		x >>= 1;
	}
	return n;
}


void createBitReversedIndices(int nn) {

	br.resize(nn);

	int log_2 = (int)log2(nn);
	for (int i = 0; i < nn; ++i) {
		br[i] = bitReverse(i, (int)log_2);
	}
}

complex w(int k, int nn, int dir) {

	float it = - dir * 2.0f * k * M_PIf / nn;
	return complex(cosf(it), sinf(it));
}


void preComputeTwiddleFactors() {

	for (int i = 0; i < FW ; ++i) {
		for (int j = 1; j <= FW ; ++j) {
			preW[i][j] = w(i, j, -1);
		}
	}
}


#define GridIn(a,b) in[(a) * FW + (b)]
#define GridOut(a,b) out[(a) * FW + (b)]

void fftStageLines(int dir, complex *in, complex *out, int iter, int log_2, int line, int column, int stage) {

	complex elemk, elemks, ww;
	int index, groupShift, shift;

	int groups = (int)(iter / pow(2.0f, stage));
	int groupSize = 2 * iter / groups;
	int k = line % (groupSize / 2);
	int group = line / (groupSize / 2);
	groupShift = (int)pow(2.0f, stage +  1);

	index = k + group * groupShift;
	shift = (int)pow(2.0f, stage);
	//printf("---- %d %d %d %d\n", index, shift, k, groupSize);

	if (stage == 0) {
		elemk = GridIn(br[index], column);
		elemks = GridIn(br[index + shift], column);
		ww = w(k, groupShift, dir) * elemks;//preW[k][groupShift]; //
		GridOut(index, column) = elemk + ww;
		GridOut(index + shift, column) = elemk - ww;
		//printf("%d %d - %f %f - %f %f - %f %f\n", index, index + shift, GridIn(index, column).r, GridIn(index + shift, column).r, GridOut(index, column).r, GridOut(index + shift, column).r, w(k, groupShift, dir).r, w(k, groupShift, dir).i);

	}
	else {
		elemk = GridIn(index, column);
		elemks = GridIn(index + shift, column);
		ww = w(k, groupShift, dir) * elemks;
		GridOut(index, column) = elemk + ww;
		GridOut(index + shift, column) = elemk - ww;
		//printf("%d %d - %f %f - %f %f - %f %f\n",index, index + shift, GridIn(index, column).r, GridIn(index + shift, column).r, GridOut(index, column).r, GridOut(index + shift, column).r, w(k, groupShift, dir).r, w(k, groupShift, dir).i);
	}
}



void fftStageColumns(int dir, complex *in, complex *out, int iter, int log_2, int line, int column, int stage) {

	complex elemk, elemks, ww;
	int index, shift;

	int p = (int)pow(2, stage);
	int groups = (int)(iter / p);
	int groupSize = 2 * iter / groups;
	int k = column % (groupSize / 2);
	int group = column / (groupSize / 2);
	int groupShift = p * 2;

	index = k + group * groupShift;
	shift = p;

	if (stage == 0) {
		elemk = GridIn(line, br[index]);
		elemks = GridIn(line, br[index + shift]);

	}
	else {
		elemk = GridIn(line, index);
		elemks = GridIn(line, index + shift);
	}
	ww = preW[k][groupShift] * elemks;
	//ww = w(k, groupShift, dir) * elemks;
	GridOut(line, index) = elemk + (ww);
	GridOut(line, index + shift) = elemk - (ww);
}


// dir = 1 => fft; dir = -1 => fftInv
complex *fft2D(int dir, complex *imgmyComplex, int width, int height) {

	complex *in1 = (complex *)malloc(sizeof(complex) * FW * FW);
	createBitReversedIndices(width);

	complex *in0 = imgmyComplex;
	complex *res = in0;

	int log_2 = (int)log2(width);

	//displayArray(in0, width, height);

	// rows
	for (int stage = 0; stage < log_2; ++stage) {
		//int j = 0;
		for (int i = 0; i < width; ++i) {
			for (int j = 0; j < height / 2; ++j) {
				if (stage % 2 == 0) {
					fftStageColumns(dir, in0, in1, width / 2, log_2, i, j, stage);
				}
				else {
					fftStageColumns(dir, in1, in0, width / 2, log_2, i, j, stage);
				}
			}
		}
		//if (stage % 2)
		//	displayArray(in0, width, height);
		//else
		//	displayArray(in1, width, height);
		//int x = 0;

		//for (int i = 0; i < 8; ++i) {
		//	if (stage % 2 == 0)
		//		printf("#%f - %f    ", in1[i*FW].r, in1[i*FW].i);
		//	else
		//		printf("#%f - %f    ", in0[i*FW].r, in0[i*FW].i);
		//}
		//printf("\n\n");
	}


	if (log_2 % 2 == 1) {
		complex *aux = in0;
		in0 = in1;
		in1 = aux;
	}

	//displayArray(in0, width, height);
	//displayArray(in1, width, height);

	// columns
	for (int stage = 0; stage < log_2; ++stage) {
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width / 2; ++i) {
				if (stage % 2 == 0)
					fftStageLines(dir, in0, in1, width / 2, log_2, i, j, stage);
				else
					fftStageLines(dir, in1, in0, width / 2, log_2, i, j, stage);
			}
		}
	}

	//displayArray(in0, width, height);
	//free(in1);
	in0 = res;
	//displayArray(in0, width, height);
	//if (dir == 1) {
	//	for (int i = 0; i < width; ++i) {
	//		for (int j = 0; j < width; ++j) {
	//			in0[i*width + j].real /= width*width;
	//			in0[i*width + j].img /= width*width;
	//		}
	//	}
	//}
	return in0;
}

// ------------------------------------------------------------
//
// GPU Stuff
//

GLuint texture0, texture1, texH0 = 0;
float ht0[FW][FW][4];
// first dimension is the layer
float textureData[4][FW][FW][4];
float texturePingPong[4][FW][FW][4];

float textureRead[4][FW][FW][4];



void sendDataToGPU() {

	for (int i = 0; i < FW; ++i) {
		for (int j = 0; j < FW; ++j) {
				textureData[0][i][j][0] = htk[i*FW + j].r;
				textureData[0][i][j][1] = htk[i*FW + j].i;
				textureData[0][i][j][2] = 0.0f;
				textureData[0][i][j][3] = 0.0f;
		}
	}

	for (int i = 0; i < FW; ++i) {
		for (int j = 0; j < FW; ++j) {
			textureData[1][i][j][0] = slopeX[i*FW + j].r;
			textureData[1][i][j][1] = slopeX[i*FW + j].i;
			textureData[1][i][j][2] = slopeZ[i*FW + j].r;
			textureData[1][i][j][3] = slopeZ[i*FW + j].i;
		}
	}
	for (int i = 0; i < FW; ++i) {
		for (int j = 0; j < FW; ++j) {
			textureData[2][i][j][0] = dX[i*FW + j].r;
			textureData[2][i][j][1] = dX[i*FW + j].i;
			textureData[2][i][j][2] = dZ[i*FW + j].r;
			textureData[2][i][j][3] = dZ[i*FW + j].i;
		}
	}

	if (!texture0) {
		glGenTextures(1, &texture0);
		glBindTexture(GL_TEXTURE_2D_ARRAY, texture0);
		glTexStorage3D(GL_TEXTURE_2D_ARRAY, 1, GL_RGBA32F, FW, FW, 4);

		glGenTextures(1, &texture1);
		glBindTexture(GL_TEXTURE_2D_ARRAY, texture1);
		glTexStorage3D(GL_TEXTURE_2D_ARRAY, 1, GL_RGBA32F, FW, FW, 4);
	}

	glBindTexture(GL_TEXTURE_2D_ARRAY, texture0);
	glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0,
						0, 0, 0, 
						FW, FW, 4, 
						GL_RGBA, GL_FLOAT, &textureData[0][0][0][0]);
	glBindTexture(GL_TEXTURE_2D_ARRAY, 0);
}


void getDataFromGPU(int t) {


	int size;

	glBindTexture(GL_TEXTURE_2D_ARRAY, t);
	glGetTexImage(GL_TEXTURE_2D_ARRAY, 0,  GL_RGBA, GL_FLOAT, &textureRead[0][0][0][0]);

	for (int i = 0; i < FW; ++i) {
		for (int j = 0; j < FW; ++j) {
			outH[i*FW + j].r = textureRead[0][i][j][0];
		}
	}

	for (int i = 0; i < FW; ++i) {
		for (int j = 0; j < FW; ++j) {
			outSX[i*FW + j].r = textureRead[1][i][j][0];
			outSZ[i*FW + j].r = textureRead[1][i][j][2];
		}
	}
	for (int i = 0; i < FW; ++i) {
		for (int j = 0; j < FW; ++j) {
			outDX[i*FW + j].r = textureRead[2][i][j][0];
			outDZ[i*FW + j].r = textureRead[2][i][j][2];
		}
	}
}


void executeShaders() {

	int k = (int)log2f(FW);


	programFFT_H.setUniform("log_width", k);
	programFFT_H.setUniform("pingpong0", 0);
	programFFT_H.setUniform("pingpong1", 1);

	glBindImageTexture(0, texture0, 0,
		GL_TRUE, 0, GL_READ_WRITE, GL_RGBA32F);
	glBindImageTexture(1, texture1, 0,
		GL_TRUE, 0, GL_READ_WRITE, GL_RGBA32F);

	glUseProgram(programFFT_H.getProgramIndex());

	for (int i = 0; i < k; ++i) {

		programFFT_H.setUniform("pingpong", i % 2);
		programFFT_H.setUniform("stage", i);

		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
		glDispatchCompute(FW / 8, FW / 8, 1);

		//glFinish();
	}
/*	if (true)
		getDataFromGPU(texture0);
	else
		getDataFromGPU(texture1);
	int x = 0;
*/
	programFFT_V.setUniform("log_width", k);
	programFFT_V.setUniform("pingpong0", 0);
	programFFT_V.setUniform("pingpong1", 1);
	glBindImageTexture(0, texture0, 0,
		GL_TRUE, 0, GL_READ_WRITE, GL_RGBA32F);
	glBindImageTexture(1, texture1, 0,
		GL_TRUE, 0, GL_READ_WRITE, GL_RGBA32F);

	glUseProgram(programFFT_V.getProgramIndex());
	for (int i = 0; i < k; ++i) {

		programFFT_V.setUniform("pingpong", (k+i) % 2);
		programFFT_V.setUniform("stage", i);

		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
		glDispatchCompute(FW / 8, FW / 8, 1);
	//	glFinish();
	}

/*	getDataFromGPU(texture0);
	x = 2;
*/
}


void sendH02GPU() {

	float *t = (float *)malloc(FW * FW * 4 * sizeof(float));

	for (int i = 0; i < FW; ++i) {
		for (int j = 0; j < FW; ++j) {
			int index = j * FW + i;
			t[index*4 + 0] = h0k[i*FW+j].r;
			t[index*4 + 1] = h0k[i*FW + j].i;
			t[index*4 + 2] = h0kconj[i*FW + j].r;
			t[index*4 + 3] = h0kconj[i*FW + j].i;
		}
	}
	if (!texH0) {
		glGenTextures(1, &texH0);
		glBindTexture(GL_TEXTURE_2D, texH0);
		glTexStorage3D(GL_TEXTURE_2D_ARRAY, 1, GL_RGBA32F, FW, FW, 1);
	}
	glBindTexture(GL_TEXTURE_2D, texH0);
	glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0,
		0, 0, 0,
		FW, FW, 1,
		GL_RGBA, GL_FLOAT, &t[0]);

	free(t);
}


// ------------------------------------------------------------
//
// Init Wave Field
//


complex getNormalDistSample() {

	float x1 = rand() * 1.0f / RAND_MAX ;
	float x2 = rand() * 1.0f / RAND_MAX ;

	x1 == 0 ? x1 = 0.00000001f: x1 = x1;
	x2 == 0 ? x2 = 0.00000001f : x2 = x2;

	float z1 = sqrtf(-2 * logf(x1))*cosf(2 * M_PIf * x2);
	float z2 = sqrtf(-2 * logf(x1))*sinf(2 * M_PIf * x2);

	return complex(z1, z2);
}


float Phillips(point2d vecK) {

	if (vecK.length() == 0.0f) 
		return 0.0f;

	float k = vecK.length();
	float l = windSpeed * windSpeed / 9.8f;
	point2d khat = vecK.normalize();

	float dot_k_w = khat.dot(windDir.normalize());
	float result = A * exp(-1 / (k*l*k*l)) / pow(k, 4) * pow(dot_k_w, 2);
	// different from shaders
	result *= expf(-k*k*0.1f*0.1f);  
	float ac = pow(2 * M_PIf / L, 2);

	return result*ac ;
}

float dispersion(float k) {

	return sqrtf(9.8f * k);
}


void initWaveField() {

	// build h0k and h0kconj
	srand(0);

	int index;

	for (int i = 0; i < FW; ++i) {
		for (int j = 0; j < FW; ++j) {

			index = i * FW + j;
			// should subtract 0.5 to get symmetry?
			point2d k = point2d(2 * M_PIf * (i - FW / 2.0f) / L, 2 * M_PIf * (j - FW / 2.0f) / L);
			complex c1 = getNormalDistSample();
			complex c2 = getNormalDistSample();

			// shouldn´t the conjugate be stored in the same array?
			h0k[index] = c1 * sqrtf(Phillips(k) * 0.5f);
			h0kconj[index] = conjugate(c2 * sqrtf(Phillips(-k) * 0.5f));
			//h0kconj[index] = conjugate(h0k[index]);

		}
	}
	sendH02GPU();

	//displayArray(h0k, FW, FW);
	//displayArray(h0kconj, FW, FW);
	outH = (complex *)fftwf_malloc(sizeof(complex) * FW * FW);
	outSX = (complex *)fftwf_malloc(sizeof(complex) * FW * FW);
	outSZ = (complex *)fftwf_malloc(sizeof(complex) * FW * FW);
	outDX = (complex *)fftwf_malloc(sizeof(complex) * FW * FW);
	outDZ = (complex *)fftwf_malloc(sizeof(complex) * FW * FW);

	preComputeTwiddleFactors();
}


void iterateWaveField() {

	int index;
	int sign;

	if (htkCPU) {
		for (int i = 0; i < FW; ++i) {
			for (int j = 0; j < FW; ++j) {

				index = i * FW + j;
				point2d k = point2d(2 * M_PIf * (i - FW / 2.0f) / L, 2 * M_PIf * (j - FW / 2.0f) / L);

				point2d knorm = k.normalize();
				float klen = k.length();
				complex t = h0k[index] * euler(dispersion(klen)*timer);
				complex tc = h0kconj[index] * euler(-dispersion(klen)*timer);

				htk[index] = (t + tc);
				slopeX[index] = complex(0, k.x) * htk[index];
				slopeZ[index] = complex(0, k.y) * htk[index];
				dX[index] = complex(0, -knorm.x) * htk[index];;
				dZ[index] = complex(0, -knorm.y) * htk[index];;
			}
		}
	}
	else {
		glUseProgram(programHTK.getProgramIndex());
		glBindImageTexture(0, texture0, 0,
			GL_TRUE, 0, GL_READ_WRITE, GL_RGBA32F);
		programHTK.setUniform("tilde_h0k", (int)texH0);
		programHTK.setUniform("tilde_hkt", 0);
		programHTK.setUniform("timer", timer);

		glDispatchCompute(FW / 8, FW / 8, 1);
	}


	if (useFFTW == 0) {

		fftwf_plan height = fftwf_plan_dft_2d(FW, FW, (fftwf_complex *)(&(htk[0])), (fftwf_complex *)&outH[0], FFTW_BACKWARD, FFTW_ESTIMATE);
		fftwf_plan sx = fftwf_plan_dft_2d(FW, FW, (fftwf_complex *)(&(slopeX[0])), (fftwf_complex *)&outSX[0], FFTW_BACKWARD, FFTW_ESTIMATE);
		fftwf_plan sz = fftwf_plan_dft_2d(FW, FW, (fftwf_complex *)(&(slopeZ[0])), (fftwf_complex *)&outSZ[0], FFTW_BACKWARD, FFTW_ESTIMATE);
		fftwf_plan dx = fftwf_plan_dft_2d(FW, FW, (fftwf_complex *)(&(dX[0])), (fftwf_complex *)&outDX[0], FFTW_BACKWARD, FFTW_ESTIMATE);
		fftwf_plan dz = fftwf_plan_dft_2d(FW, FW, (fftwf_complex *)(&(dZ[0])), (fftwf_complex *)&outDZ[0], FFTW_BACKWARD, FFTW_ESTIMATE);
		fftwf_execute(height);
		fftwf_execute(sx);
		fftwf_execute(sz);
		fftwf_execute(dx);
		fftwf_execute(dz);
		fftwf_destroy_plan(height);
		fftwf_destroy_plan(sx);
		fftwf_destroy_plan(sz);
		fftwf_destroy_plan(dx);
		fftwf_destroy_plan(dz);
	}
	else if (useFFTW == 1) {
		FFT2D(htk, FW, FW, -1);
		outH = htk;
		FFT2D(slopeX, FW, FW, -1);
		outSX = slopeX;
		FFT2D(slopeZ, FW, FW, -1);
		outSZ = slopeZ;
		FFT2D(dX, FW, FW, -1);
		outDX = dX;
		FFT2D(dZ, FW, FW, -1);
		outDZ = dZ;
	}
	else if (useFFTW == 2) {
		outH = fft2D(-1, htk, FW, FW);
		outSX = fft2D(-1, slopeX, FW, FW);
		outSZ = fft2D(-1, slopeZ, FW, FW);
		outDX = fft2D(-1, dX, FW, FW);
		outDZ = fft2D(-1, dZ, FW, FW);
	}
	else {
		if (htkCPU)
			sendDataToGPU();
		{
			PROFILE_GL("shaders");
			executeShaders();
		}
		getDataFromGPU(texture0);
	}


	for (int i = 0; i < FW; ++i) {
		for (int j = 0; j < FW; ++j) {

			index = i * FW + j;
			if ((i + j) % 2)
				sign = -1;
			else sign = 1;

			float scale = 1;

			pos[index * 4 + 0] = (i - FW / 2.0f) * L  * 1.0f / FW - sign * outDX[index].r ;
			pos[index * 4 + 1] = sign * outH[index].r;
			pos[index * 4 + 2] = (j - FW / 2.0f) * L  * 1.0f / FW - sign * outDZ[index].r;
			pos[index * 4 + 3] = 1.0f;

			normal[index * 3 + 0] = -sign * outSX[index].r / scale;
			normal[index * 3 + 1] = 1;
			normal[index * 3 + 2] = -sign * outSZ[index].r / scale;

			if (sign * outH[index].r > max)
				max = sign * outH[index].r;
			if (sign * outH[index].r <min)
				min = sign * outH[index].r;
		}
	}



	int bufferPos = sg.mMyMeshes[0].vboPos;
	glBindBuffer(GL_ARRAY_BUFFER, bufferPos);
	glBufferData(GL_ARRAY_BUFFER, FW * FW * 4 * sizeof(float), pos, GL_DYNAMIC_DRAW);

	int bufferNormal = sg.mMyMeshes[0].vboNormal;
	glBindBuffer(GL_ARRAY_BUFFER, bufferNormal);
	glBufferData(GL_ARRAY_BUFFER, FW * FW * 3 * sizeof(float), normal, GL_DYNAMIC_DRAW);
}

// ------------------------------------------------------------
//
// Reshape Callback Function
//

void changeSize(int w, int h) {

	float ratio;
	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if(h == 0)
		h = 1;

	// Set the viewport to be the entire window
    glViewport(0, 0, w, h);

	ratio = (1.0f * w) / h;
	vsml->loadIdentity(VSMathLib::PROJECTION);
	vsml->perspective(53.13f, ratio, 0.1f, 10000.0f);
	//vsml->ortho(-2 , 2 , -2/ratio, 2/ratio, -10, 10);
}


// ------------------------------------------------------------
//
// Render stuff
//
void renderScene(void) {

	{
		PROFILE("Frame");
		{
			PROFILE("fft");
			iterateWaveField();
		}
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Set both matrices to the identity matrix
		vsml->loadIdentity(VSMathLib::VIEW);
		vsml->loadIdentity(VSMathLib::MODEL);

		// set camera
		vsml->lookAt(camX, camY, camZ, 0,0,0, 0,1,0);

		float res[4];
		vsml->multMatrixPoint(VSMathLib::VIEW, lightDir, res);
		vsml->normalize(res);
		program.setBlockUniform("Lights", "l_dir", res);

		{
			PROFILE_GL("Render models");

			// set the shader to render models
			glUseProgram(program.getProgramIndex());
			// start counting primitives
			glBeginQuery(GL_PRIMITIVES_GENERATED, counterQ);
			axis.render();

			//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			sg.render();
			//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);


			// stop counting primitives
			glEndQuery(GL_PRIMITIVES_GENERATED);
		}

		// FPS computation and display
		frame++;
		myTime=glutGet(GLUT_ELAPSED_TIME);
		if (myTime - timebase > 1000) {
				sprintf(s,"FPS:%4.2f  Triangles: %d  L: %d Max: %f   Min %f",
					frame*1000.0/(myTime-timebase) , primitiveCounter, L, max, min);
			timebase = myTime;
			frame = 0;
#if (__VSL_FONT_LOADING__ == 1)
			vsfl.prepareSentence(aSentence,s);
#endif
		}

#if (__VSL_FONT_LOADING__ == 1)
		// Display text info
		{
			PROFILE("Dump");
			//set the shader for rendering the sentence
			glUseProgram(programFonts.getProgramIndex());
			// prepare sentence with profile info
			std::string s = VSProfileLib::DumpLevels();
			vsfl.prepareSentence(profileSentence, s);
			//set the shader for rendering the sentence
			// render sentences
			vsfl.renderSentence(10,10,aSentence);
			vsfl.renderSentence(10, 30, profileSentence);

		}
#endif
		 //swap buffers
		{
			PROFILE("Swap");
			glutSwapBuffers();
		}
	} // end PROFILE("Frame")
	{
		PROFILE("Collect GL Queries Time");
		VSProfileLib::CollectQueryResults();
		glGetQueryObjectuiv(counterQ, GL_QUERY_RESULT, &primitiveCounter);
	}

	timer += 0.016f;
}


// ------------------------------------------------------------
//
// Events from the Keyboard
//

void processKeys(unsigned char key, int xx, int yy)
{
	std::string s;
	switch(key) {

		case 27:

			glutLeaveMainLoop();
			break;

		case 'z': r -= 0.1f;
				camX = r * sin(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
				camZ = r * cos(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
				camY = r *   						     sin(beta * 3.14f / 180.0f);
				break;
		case 'x': r += 0.1f;
				camX = r * sin(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
				camZ = r * cos(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
				camY = r *   						     sin(beta * 3.14f / 180.0f);
				break;
		case 'm': glEnable(GL_MULTISAMPLE); break;
		case 'n': glDisable(GL_MULTISAMPLE); break;
		case 'k': VSProfileLib::Reset(); break;
		case 'p': s = VSProfileLib::DumpLevels();
				printf("%s\n", s.c_str());
				break;
		case '+': L *= 2; max = -1000; min = 1000;
			programHTK.setUniform("L", L);
			initWaveField();
			break;
		case '-': L /= 2; max = -1000; min = 1000;
			programHTK.setUniform("L", L);
			initWaveField();
			break;
		case 'f': useFFTW = (++useFFTW) % 4;
			//for (int i = 0; i < 10; ++i) {
			//	printf("%f %f\n", htk[i].r, htk[i].i);
			//}
			printf("\nFFT: %d \n", useFFTW);
			VSProfileLib::Reset();
			//initWaveField();
			//timer = 0;
			break;
		case 'h': htkCPU = !htkCPU;
			printf("\nhtkCPU: %d \n", htkCPU);
			VSProfileLib::Reset();
			break;

	}
	camX = r * sin(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
	camZ = r * cos(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
	camY = r *   						     sin(beta * 3.14f / 180.0f);

//  uncomment this if not using an idle func
//	glutPostRedisplay();
}


// ------------------------------------------------------------
//
// Mouse Events
//

void processMouseButtons(int button, int state, int xx, int yy)
{
	// start tracking the mouse
	if (state == GLUT_DOWN)  {
		startX = xx;
		startY = yy;
		if (button == GLUT_LEFT_BUTTON)
			tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON)
			tracking = 2;
	}

	//stop tracking the mouse
	else if (state == GLUT_UP) {
		if (tracking == 1) {
			alpha -= (xx - startX);
			beta += (yy - startY);
		}
		else if (tracking == 2) {
			r += (yy - startY) * 0.1f;
		}
		tracking = 0;
	}
}

// Track mouse motion while buttons are pressed
void processMouseMotion(int xx, int yy)
{

	int deltaX, deltaY;
	float alphaAux, betaAux;
	float rAux;

	deltaX =  - xx + startX;
	deltaY =    yy - startY;

	// left mouse button: move camera
	if (tracking == 1) {


		alphaAux = alpha + deltaX;
		betaAux = beta + deltaY;

		if (betaAux > 85.0f)
			betaAux = 85.0f;
		else if (betaAux < -85.0f)
			betaAux = -85.0f;

		rAux = r;

		camX = rAux * sin(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
		camZ = rAux * cos(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
		camY = rAux *   						       sin(betaAux * 3.14f / 180.0f);
	}
	// right mouse button: zoom
	else if (tracking == 2) {

		alphaAux = alpha;
		betaAux = beta;
		rAux = r + (deltaY * 0.1f);

		camX = rAux * sin(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
		camZ = rAux * cos(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
		camY = rAux *   						       sin(betaAux * 3.14f / 180.0f);
	}


//  uncomment this if not using an idle func
//	glutPostRedisplay();
}


void mouseWheel(int wheel, int direction, int x, int y) {

	r += direction * 0.1f;
	camX = r * sin(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
	camZ = r * cos(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
	camY = r *   						     sin(beta * 3.14f / 180.0f);

//  uncomment this if not using an idle func
//	glutPostRedisplay();
}


// --------------------------------------------------------
//
// Shader Stuff
//


GLuint setupShaders() {

	std::string path = PATH_TO_FILES;

	// shader for FFTs
	std::vector<std::string> fnH = { path + "shaders/header.glsl", path + "shaders/defines.glsl", path + "shaders/definesGlobal.glsl", path + "shaders/butterflyH.comp"};
	programFFT_H.init();
	programFFT_H.loadShader(VSShaderLib::COMPUTE_SHADER, fnH);

	programFFT_H.prepareProgram();

	VSGLInfoLib::getUniformsInfo(programFFT_H.getProgramIndex());

	programFFT_H.setUniform("log_width", (int)log2(FW));

	printf("InfoLog for FFT_H Shader\n%s\n\n", programFFT_H.getAllInfoLogs().c_str());

	std::vector<std::string> fnV = { path + "shaders/header.glsl", path + "shaders/defines.glsl", path + "shaders/definesGlobal.glsl", path + "shaders/butterflyV.comp" };
	programFFT_V.init();
	programFFT_V.loadShader(VSShaderLib::COMPUTE_SHADER, fnV);

	programFFT_V.prepareProgram();

	programFFT_V.setUniform("log_width", (int)log2(FW));

	VSGLInfoLib::getUniformsInfo(programFFT_V.getProgramIndex());
	printf("InfoLog for FFT_V Shader\n%s\n\n", programFFT_V.getAllInfoLogs().c_str());


	// Shader for htk
	std::vector<std::string> fnHTK = { path + "shaders/header.glsl", path + "shaders/defines.glsl", path + "shaders/definesGlobal.glsl", 
									 path + "shaders/directional.glsl", path + "shaders/fft_hkt.comp" };

	programHTK.init();
	programHTK.loadShader(VSShaderLib::COMPUTE_SHADER, fnHTK);
	programHTK.prepareProgram();

	programHTK.setUniform("width", FW);
	programHTK.setUniform("dispersionMode", 0);
	programHTK.setUniform("depth", 100);




	VSGLInfoLib::getUniformsInfo(programHTK.getProgramIndex());
	printf("InfoLog for HTK Shader\n%s\n\n", programHTK.getAllInfoLogs().c_str());


	// Shader for fonts
	programFonts.init();
	programFonts.loadShader(VSShaderLib::VERTEX_SHADER, path + "shaders/color.vert");
	programFonts.loadShader(VSShaderLib::FRAGMENT_SHADER, path + "shaders/color.frag");

	// set semantics for the shader variables
	programFonts.setProgramOutput(0,"outputF");
	programFonts.setVertexAttribName(VSShaderLib::VERTEX_COORD_ATTRIB, "position");
	programFonts.setVertexAttribName(VSShaderLib::TEXTURE_COORD_ATTRIB, "texCoord");


	programFonts.prepareProgram();
	VSGLInfoLib::getUniformsInfo(programFonts.getProgramIndex());

	// add sampler uniforms
	programFonts.setUniform("texUnit", 0);

	printf("InfoLog for Font Shader\n%s\n\n", programFonts.getAllInfoLogs().c_str());


	// Shader formodels
	program.init();
	program.loadShader(VSShaderLib::VERTEX_SHADER, path + "shaders/pixeldirdifambspec.vert");
	program.loadShader(VSShaderLib::FRAGMENT_SHADER, path + "shaders/pixeldirdifambspec.frag");

	// set semantics for the shader variables
	program.setProgramOutput(0, "colorOut");
	program.setVertexAttribName(VSShaderLib::VERTEX_COORD_ATTRIB, "position");
	program.setVertexAttribName(VSShaderLib::TEXTURE_COORD_ATTRIB, "texCoord");
	program.setVertexAttribName(VSShaderLib::NORMAL_ATTRIB, "normal");

	program.prepareProgram();

	VSGLInfoLib::getProgramInfo(program.getProgramIndex());
	VSGLInfoLib::getUniformsInfo(program.getProgramIndex());

	printf("InfoLog for Model Shader\n%s\n", program.getAllInfoLogs().c_str());
	// set sampler uniform
	program.setUniform("texUnit", 0);

	return program.isProgramValid();
}



// ------------------------------------------------------------
//
// Model loading and OpenGL setup
//


int init() {

	axis.set(5, 0.02f);

	float blue[4] = { 0.3f, 0.4f, 0.8f, 1.0f };
	float white[4] = { 0.8f, 0.8f, 0.8f, 1.0f };
	float shininess = 64;
	sg.set(FW * 0.5f, FW);
	sg.setColor(VSResourceLib::DIFFUSE, blue);
	sg.setColor(VSResourceLib::SPECULAR, white);
	sg.setColor(VSResourceLib::SHININESS, &shininess);

		
	gridY.set(VSGrid::Y, 5, 25);
	// some GL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_MULTISAMPLE);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClearColor(0.75f, 0.75f, 0.75f, 0.75f);
	glClearColor(0.25f, 0.25f, 0.25f, 0.25f);
	// generate a query to count primitives
	glGenQueries(1,&counterQ);

	return true;

}


void initVSL() {

	// set the material's block name
	VSResourceLib::setMaterialBlockName("Material");

	// Init VSML
	vsml = VSMathLib::getInstance();
	vsml->setUniformBlockName("Matrices");	
	vsml->setUniformName(VSMathLib::PROJ_VIEW_MODEL, "m_pvm");
	vsml->setUniformName(VSMathLib::NORMAL, "m_normal");
	vsml->setUniformName(VSMathLib::VIEW_MODEL, "m_viewModel");


#if (__VSL_TEXTURE_LOADING__ == 1)

	// Init VSFL Fonts
	std::string path = PATH_TO_FILES;
	vsfl.load(path + "fonts/couriernew10");
	vsfl.setFixedFont(true);
	vsfl.setColor(1.0f, 0.5f, 0.25f, 1.0f);
	aSentence = vsfl.genSentence();
	profileSentence = vsfl.genSentence();
#endif
}



// ------------------------------------------------------------
//
// Main function
//


int main(int argc, char **argv) {

//  GLUT initialization
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA|GLUT_MULTISAMPLE);

	// Set context
	glutInitContextVersion (4, 6);
	glutInitContextProfile (GLUT_CORE_PROFILE );

	glutInitWindowPosition(100,100);
	glutInitWindowSize(640,360);
	glutCreateWindow("Lighthouse3D - VSL Demo");


//  Callback Registration
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutIdleFunc(renderScene);

//	Mouse and Keyboard Callbacks
	glutKeyboardFunc(processKeys);
	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);
	glutMouseWheelFunc ( mouseWheel ) ;

//	return from main loop
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

//	Init GLEW
	glewExperimental = GL_TRUE;
	glewInit();
	if (!glewIsSupported("GL_VERSION_4_5")) {
		printf("OpenGL 4.6 not supported\n");
		exit(1);
	}

	// this function will display generic OpenGL information
	VSGLInfoLib::getGeneralInfo();

	initVSL();


	// init OpenGL 
	if (!init()) {
		printf("Some init problem!!!\n");
		exit(1);
	}

	setupShaders();
	initWaveField();

	//iterateWaveField();
	//for (int i = 0; i < 10; ++i) {
	//	printf("%f %f\n", outH[i].r, outH[i].i);
	//}
	//printf("\n");

	//useFFTW = false;

	//initWaveField();

	//iterateWaveField();
	//for (int i = 0; i < 10; ++i) {
	//	printf("%f %f\n", outH[i].r, outH[i].i);
	//}
	//printf("\n");

	//sendDataToGPU();
	//getDataFromGPU();


	//  GLUT main loop
	glutMainLoop();

	return(1);

}

