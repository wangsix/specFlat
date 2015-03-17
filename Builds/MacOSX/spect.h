//
//  spect.h
//  SpecFlat
//
//  Created by Cheng-i Wang on 12/3/13.
//
//	FFT codes adapted from Tom Erbe`s code from MUS267 class FALL2013@UCSD

#ifndef __SpecFlat__spect__
#define __SpecFlat__spect__
#if WINDOWS
#define FFT_USED_FFTW 1
#define FFT_USED_APPLEVECLIB 0
#define FFT_USED_INTEL 0
#else
#define FFT_USED_FFTW 0
#define FFT_USED_APPLEVECLIB 1
#define FFT_USED_INTEL 0
#endif


#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if FFT_USED_FFTW
#include <rfftw.h>
#elif FFT_USED_APPLEVECLIB
#include <Accelerate/Accelerate.h>
#elif FFT_USED_INTEL
#endif

#define TEndian16_Swap(value)                 \
(((((unsigned short)value)<<8) & 0xFF00)   | \
((((unsigned short)value)>>8) & 0x00FF))

#define TEndian32_Swap(value)                     \
(((((unsigned long)value)<<24) & 0xFF000000)  | \
((((unsigned long)value)<< 8) & 0x00FF0000)  | \
((((unsigned long)value)>> 8) & 0x0000FF00)  | \
((((unsigned long)value)>>24) & 0x000000FF))

struct swap_4
{
	unsigned char b1, b2, b3, b4;
};

static inline float swapfloat (float f)
{
	unsigned char tmp;
	swap_4 *pf = (swap_4 *)&f;
	tmp = pf->b1; pf->b1 = pf->b4; pf->b4 = tmp;
	tmp = pf->b2; pf->b2 = pf->b3; pf->b3 = tmp;
	return f;
}

enum
{
	kSizeFFT = 1024,
	kHalfSizeFFT = 512,
	kMaxSizeFFT = 65536,
	kHalfMaxSizeFFT = 32768
};
enum
{
	kMonoMode,
	kStereoMode,
	kMono2StereoMode
};
enum
{
	kHamming,
	kVonHann,
	kBlackman,
	kKaiser
};

class spectFlat{
	
public:
	spectFlat();
	~spectFlat();
	
	bool	AllocateMemory();
	void	FreeMemory();
	void	Process(float **data, long sampleFrames);
	void	Suspend();
	void	SetBypass(float b);
	void	SetInNumChannels(long n);
	void	SetOutNumChannels(long n);
	void	SetMinOrMax(float m);
	void	SetWhiteOrPink(float m);
	void	SetStepSize(float size);
	void	SetChannelMode(long nIn, long nOut);
	void	SetDecaySize(float n);
	bool	GetMinOrMax(void){return max;}; // Did not use.
	float	GetStepSize(void){return stepSize;}; // Did not use.
	float	*GetWeights(void){return weights;}; // Did not use.
	void	initHammingWindows(void);
	void	initVonHannWindows(void); // Did not implement.
	void	initBlackmanWindows(void); // Did not implement.
	void	initKaiserWindows(void); // Did not implement.
	void	scaleWindows(void);
		
private:
	void	ProcessBlockActual();
	void	ProcessSpect();
	
	// The manipulation is implemented here: --------------------------------
	void	cartToPolar(float *spectrum, float *polarSpectrum, float *decay);
	void	polarToCart(float *polarSpectrum, float *spectrum, float *weight, float *decay);
	//-----------------------------------------------------------------------
	
	long	channelMode;
	long	winType;
	long	nMaxChannels;
	long	nInChannels, nOutChannels;
	long	halfSizeFFT, sizeFFT, log2n, logScale;
	long	blockSize, bufferPosition;
	long	inputTimeL, outputTimeL;
	long	inputTimeR, outputTimeR;
	float	oneOverBlockSize;
	float	oneOverFFTSize;
	
	float	bypass;
	float	max;
	float	white;
	float	stepSize;
	float	decaySize;
	
	float	pi, twoPi;
	float	magSum;
	float	*weights;
	float	**weightsPt;
	float	*analysisWindow, *synthesisWindow;
	float	*inBuffer, *outBuffer, *inFFT, *polarSpect, *decaySpect;
	float	*inShift, *outShift, *inSpectral, *outSpectral, *outDisplaySpect;
	float	**inBufferPt, **outBufferPt, **inFFTPt, **polarSpectPt, **decaySpectPt;
	float	**inShiftPt, **outShiftPt, **inSpectralPt, **outSpectralPt, **outDisplaySpectPt;
	
#if FFT_USED_FFTW
	rfftw_plan planIn, planOut;
#elif FFT_USED_APPLEVECLIB
	FFTSetup setupReal;
	COMPLEX_SPLIT split;
#elif FFT_USED_INTEL
#endif
};

#endif





