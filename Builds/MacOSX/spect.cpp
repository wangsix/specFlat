//
//  spect.cpp
//  SpecFlat
//
//  Created by Cheng-i Wang on 12/3/13.
//
//	FFT codes adapted from Tom Erbe`s code from MUS267 class FALL2013@UCSD

#include "spect.h"

spectFlat::spectFlat(){
	
	outBuffer = inBuffer = inFFT = inShift = outShift = 0;
	inSpectral = outSpectral = outDisplaySpect = polarSpect = decaySpect = 0;
	outBufferPt = inBufferPt = inFFTPt = inShiftPt = outShiftPt = 0;
	inSpectralPt = outSpectralPt = outDisplaySpectPt = polarSpectPt = decaySpectPt = 0;
	synthesisWindow = analysisWindow = weights = 0;

	bufferPosition = 0;
	inputTimeL = outputTimeL = 0;
	inputTimeR = outputTimeR = 0;

	sizeFFT = kSizeFFT;
	blockSize = sizeFFT >> 3;
	halfSizeFFT = sizeFFT >> 1;
	oneOverBlockSize = 1.0f/(float)blockSize;
	oneOverFFTSize = 1.0f/(float)sizeFFT;
	
	nMaxChannels = 2;

	SetMinOrMax(1.0f);
	SetStepSize(0.01f);
	SetBypass(0.0f);
	
	winType = kHamming;
	nInChannels = 2;
	nOutChannels = 2;
	
	AllocateMemory();
#if FFT_USED_FFTW
	planIn = rfftw_create_plan(sizeFFT, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
	planOut = rfftw_create_plan(sizeFFT, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
#elif FFT_USED_APPLEVECLIB
	log2n = (long)log2f((float)sizeFFT);
    setupReal = vDSP_create_fftsetup(log2n, FFT_RADIX2);
#elif FFT_USED_INTEL
#endif

	pi = 4.0f * atanf(1.0f);
	twoPi = 8.0f * atanf(1.0f);
	magSum = 0.0f;
	
	logScale = false;
	initHammingWindows();
	scaleWindows();
}

spectFlat::~spectFlat(){
	FreeMemory();
#if FFT_USED_FFTW
	rfftw_destroy_plan(planIn);
	rfftw_destroy_plan(planOut);
#elif FFT_USED_APPLEVECLIB
	vDSP_destroy_fftsetup(setupReal);
#elif FFT_USED_INTEL
#endif

}

bool spectFlat::AllocateMemory(){

#if FFT_USED_APPLEVECLIB
	split.realp = 0;
	split.realp = (float *) malloc(kHalfMaxSizeFFT*sizeof(float));
	split.imagp = 0;
	split.imagp = (float *) malloc(kHalfMaxSizeFFT*sizeof(float));
#endif
	
	outBuffer = inBuffer = inFFT = inShift = outShift = 0;
	inSpectral = outSpectral = outDisplaySpect = polarSpect = decaySpect = 0;
	outBufferPt = inBufferPt = inFFTPt = inShiftPt = outShiftPt = 0;
	inSpectralPt = outSpectralPt = outDisplaySpectPt = polarSpectPt = decaySpectPt = 0;
	synthesisWindow = analysisWindow = 0;
	
	decaySpect = (float *) malloc(nMaxChannels*halfSizeFFT*sizeof(float));
	polarSpect = (float *) malloc(nMaxChannels*kSizeFFT*sizeof(float));
	inBuffer = (float *)malloc(nMaxChannels*kMaxSizeFFT*sizeof(float));
	outBuffer = (float *)malloc(nMaxChannels*kMaxSizeFFT*sizeof(float));
	inFFT = (float *)malloc(nMaxChannels*kMaxSizeFFT*sizeof(float));
	inShift = (float *)malloc(nMaxChannels*kMaxSizeFFT*sizeof(float));
	outShift = (float *)malloc(nMaxChannels*kMaxSizeFFT*sizeof(float));
	inSpectral = (float *)malloc(nMaxChannels*kMaxSizeFFT*sizeof(float));
	outSpectral = (float *)malloc(nMaxChannels*kMaxSizeFFT*sizeof(float));
	outDisplaySpect = (float *)malloc(nMaxChannels*kMaxSizeFFT*sizeof(float));

	memset(decaySpect, 0, nMaxChannels*halfSizeFFT*sizeof(float));
	memset(inBuffer, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(outBuffer, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(inFFT, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(inShift, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(outShift, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(inSpectral, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(outSpectral, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(outDisplaySpect, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(polarSpect, 0, nMaxChannels*kSizeFFT*sizeof(float));

	decaySpectPt = (float **)malloc(nMaxChannels*sizeof(float *));
	inBufferPt = (float **)malloc(nMaxChannels*sizeof(float *));
	outBufferPt = (float **)malloc(nMaxChannels*sizeof(float *));
	inFFTPt = (float **)malloc(nMaxChannels*sizeof(float *));
	inShiftPt = (float **)malloc(nMaxChannels*sizeof(float *));
	outShiftPt = (float **)malloc(nMaxChannels*sizeof(float *));
	inSpectralPt = (float **)malloc(nMaxChannels*sizeof(float *));
	outSpectralPt = (float **)malloc(nMaxChannels*sizeof(float *));
	outDisplaySpectPt = (float **)malloc(nMaxChannels*sizeof(float *));
	polarSpectPt = (float **)malloc(nMaxChannels*sizeof(float *));

	weights = 0;
	weights = (float *) malloc(nMaxChannels*halfSizeFFT*sizeof(float));
	for (int i = 0; i < nMaxChannels*halfSizeFFT; i++) {
		weights[i] = 1.0f;
	}
	weightsPt = (float **)malloc(nMaxChannels*sizeof(float *));

	for (int i = 0; i < nMaxChannels; i++) {
		inBufferPt[i] = inBuffer + i*kMaxSizeFFT;
		outBufferPt[i] = outBuffer + i*kMaxSizeFFT;
		inFFTPt[i] = inFFT + i*kMaxSizeFFT;
		inShiftPt[i] = inShift + i*kMaxSizeFFT;
		outShiftPt[i] = outShift + i*kMaxSizeFFT;
		inSpectralPt[i] = inSpectral + i*kMaxSizeFFT;
		outSpectralPt[i] = outSpectral + i*kMaxSizeFFT;
		outDisplaySpectPt[i] = outDisplaySpect + i*kMaxSizeFFT;
		polarSpectPt[i] = polarSpect + i*kSizeFFT;
		decaySpectPt[i] = decaySpect + i*halfSizeFFT;
		weightsPt[i] = weights + i*halfSizeFFT;
	}
	
	synthesisWindow = 0;
	synthesisWindow = (float *) malloc(kMaxSizeFFT*sizeof(float));
	analysisWindow = 0;
	analysisWindow = (float *) malloc(kMaxSizeFFT*sizeof(float));
	
	return(true);
}

void spectFlat::FreeMemory(){
	
	for (int i = 0; i < nMaxChannels; i++) {
		if (inBufferPt[i]) {free(inBufferPt[i]); inBufferPt[i]=0;}
		if (outBufferPt[i]) {free(outBufferPt[i]); outBufferPt[i]=0;}
		if (inFFTPt[i]) {free(inFFTPt[i]); inFFTPt[i]=0;}
		if (inShiftPt[i]) {free(inShiftPt[i]); inShiftPt[i]=0;}
		if (outShiftPt[i]) {free(outShiftPt[i]); outShiftPt[i]=0;}
		if (inSpectralPt[i]) {free(inSpectralPt[i]); inSpectralPt[i]=0;}
		if (outSpectralPt[i]) {free(outSpectralPt[i]); outSpectralPt[i]=0;}
		if (outDisplaySpectPt[i]) {free(outDisplaySpectPt[i]); outDisplaySpectPt[i]=0;}
		if (polarSpectPt[i]){free(polarSpectPt[i]); polarSpectPt[i]=0;}
		if (decaySpectPt[i]){free(decaySpectPt[i]); decaySpectPt[i]=0;}
	}
	
	if(inBuffer) {free(inBuffer); inBuffer = 0; }
	if(inFFT) {free(inFFT); inFFT = 0;}
	if(outBuffer) {free(outBuffer); outBuffer = 0;}
	if(inShift) {free(inShift); inShift = 0;}
	if(outShift) {free(outShift);outShift = 0;}
	if(inSpectral) {free(inSpectral);inSpectral = 0;}
	if(outSpectral) {free(outSpectral);outSpectral = 0;}
	if(outDisplaySpect) {free(outDisplaySpect);outDisplaySpect = 0;}
	if(polarSpect) {free(polarSpect);polarSpect = 0;}
	if(decaySpect) {free(decaySpect);decaySpect = 0;}
	
	if(analysisWindow) {free(analysisWindow); analysisWindow = 0;}
	if(synthesisWindow) {free(synthesisWindow); synthesisWindow = 0;}
	if(weights) {free(weights); weights = 0;}
	
#if FFT_USED_APPLEVECLIB
	if(split.realp) {free(split.realp); split.realp = 0;}
	if(split.imagp) {free(split.imagp); split.imagp = 0;}
#endif
}

void spectFlat::Suspend(){
	memset(inBuffer, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(outBuffer, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(inFFT, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(inShift, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(outShift, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(inSpectral, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(outSpectral, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(outDisplaySpect, 0, nMaxChannels*kMaxSizeFFT*sizeof(float));
	memset(polarSpect, 0, nMaxChannels*kSizeFFT*sizeof(float));
	memset(decaySpect, 0, nMaxChannels*halfSizeFFT*sizeof(float));

	for (int i = 0; i < nMaxChannels*halfSizeFFT; i++) {
		weights[i] = 1.0f;
	}
}

void spectFlat::SetInNumChannels(long n){
	
	if (n > nMaxChannels) {
		nInChannels = nMaxChannels;
	}else if (n < 1){
		nInChannels = 1;
	}else{
		nInChannels = n;
	}
}

void spectFlat::SetOutNumChannels(long n){
	
	if (n > nMaxChannels) {
		nOutChannels = nMaxChannels;
	}else if (n < 1){
		nOutChannels = 1;
	}else{
		nOutChannels = n;
	}
}

void spectFlat::SetBypass(float b){
	if (b!=0.0f) {
		bypass = 1.0f;
	}else{
		bypass = 0.0f;
	}
}

void spectFlat::SetMinOrMax(float m){
	// True for max, false for min.
	if (m > 0.0f) {max = 1.0f;}
	else{max = -1.0f;}
}

void spectFlat::SetWhiteOrPink(float m){
	if (m > 0.0f) { white = 1.0f;}
	else{white = 0.0f;}
}

void spectFlat::SetStepSize(float step){
	
	stepSize=step;
	if(stepSize > 1.0f){
		stepSize = 1.0f;
	}
	else if (stepSize <= 0.01f){
		stepSize = 0.01f;
	}
}

void spectFlat::SetChannelMode(long nIn, long nOut){

	SetInNumChannels(nIn);
	SetOutNumChannels(nOut);
	
	if(nInChannels == 2){
		channelMode = kStereoMode;
	}else if(nOutChannels == 2){
		channelMode = kMono2StereoMode;
	}else{
		channelMode = kMonoMode;		
	}
}

void spectFlat::SetDecaySize(float n){
	if (n < 1.0f) {
		decaySize = 1.0f;
	}else{
		decaySize = n;
	}
}

void spectFlat::Process(float **data, long sampleFrames){
	
	long i, framesLeft, processframes;
	float *in1 = data[0];
	float *in2 = data[1];
	float *out1 = data[0];
	float *out2 = data[1];
		
	framesLeft = sampleFrames;
	while (framesLeft > 0)
	{
		// how many frames can we process now
		// with this we insure that we stop on the
		// blockSize boundary
		if(framesLeft+bufferPosition < blockSize)
			processframes = framesLeft;
		else
			processframes = blockSize - bufferPosition;
		// flush out the previous output, copy in the new input...
		memcpy(inBufferPt[0]+bufferPosition, in1, processframes*sizeof(float));
		for(i=0; i<processframes; i++)
		{
				out1[i] = outBufferPt[0][i+bufferPosition];
		}
		if(channelMode == kStereoMode)
		{
			memcpy(inBufferPt[1]+bufferPosition, in2, processframes*sizeof(float));
			for(i=0; i<processframes; i++)
			{
				// copy the old output into the out buffers
				out2[i] = outBufferPt[1][i+bufferPosition];
			}
			in2 += processframes;
			out2 += processframes;
		}
		else if(channelMode == kMono2StereoMode)
		{
			for(i=0; i<processframes; i++)
			{
				// copy the old output into the out buffers
				out2[i] = outBufferPt[0][i+bufferPosition];
			}
			out2 += processframes;
		}
		
		bufferPosition += processframes;
		// if filled a buffer, we process a new block
		if(bufferPosition >= blockSize)
		{
			bufferPosition = 0;
			if(bypass == true)
			{
				memcpy(outBufferPt[0], inBufferPt[0], (blockSize) * sizeof(float));
				if(channelMode == kStereoMode)
					memcpy(outBufferPt[1], inBufferPt[1], (blockSize) * sizeof(float));
			}
			else
				ProcessBlockActual();
		}
		in1 += processframes;
		out1 += processframes;
		framesLeft -= processframes;
	}
}

void spectFlat::ProcessBlockActual(){
	long	i, j;
	float	outTemp;
	long	maskFFT;
	
	maskFFT = sizeFFT - 1;
	//
	inputTimeL += blockSize;
	inputTimeR += blockSize;
	outputTimeL += blockSize;
	outputTimeR += blockSize;
	inputTimeL = inputTimeL & maskFFT;
	inputTimeR = inputTimeR & maskFFT;
	outputTimeL = outputTimeL & maskFFT;
	outputTimeR = outputTimeR & maskFFT;
	//
	// a - shift output buffer and zero out new location
	memcpy(outBufferPt[0], outBufferPt[0]+blockSize, (sizeFFT - blockSize) * sizeof(float));
	memset(outBufferPt[0]+(sizeFFT - blockSize), 0, blockSize * sizeof(float));
	if(channelMode == kStereoMode)
	{
		// a - shift output buffer and zero out new location
		memcpy(outBufferPt[1], outBufferPt[1]+blockSize, (sizeFFT - blockSize) * sizeof(float));
		memset(outBufferPt[1]+(sizeFFT - blockSize), 0, blockSize * sizeof(float));
	}
	// a - shift current samples in inShift
	memcpy(inShiftPt[0], inShiftPt[0]+blockSize, (sizeFFT - blockSize) * sizeof(float));
	// a1 - copy the new stuff in
	memcpy(inShiftPt[0]+(sizeFFT - blockSize), inBufferPt[0], blockSize * sizeof(float));
	// b - window the block in
	for(i = 0; i < sizeFFT; i++)
	{
		*(inFFTPt[0] + inputTimeL) = *(inShiftPt[0] + i) * *(analysisWindow + i);
		++inputTimeL;
		inputTimeL = inputTimeL & maskFFT;
	}
#if FFT_USED_FFTW
	rfftw_one(planIn, inFFTPt[0], inSpectralPt[0]);
#elif FFT_USED_APPLEVECLIB
	vDSP_ctoz ((COMPLEX *)inFFTPt[0], 2, &split, 1, halfSizeFFT);
	vDSP_fft_zrip (setupReal, &split, 1, log2n, FFT_FORWARD);
	memcpy(inSpectralPt[0], split.realp, halfSizeFFT * sizeof(float));
	*(inSpectralPt[0] + halfSizeFFT) = *(split.imagp);
	for(i = halfSizeFFT + 1, j = halfSizeFFT - 1; i < sizeFFT; i++, j--)
		*(inSpectralPt[0] + i) = *(split.imagp + j);
#elif FFT_USED_INTEL
#endif
	if(channelMode == kStereoMode)
	{
		// a - shift current samples in inShift
		memcpy(inShiftPt[1], inShiftPt[1]+blockSize, (sizeFFT - blockSize) * sizeof(float));
		// a1 - copy the new stuff in
		memcpy(inShiftPt[1]+(sizeFFT - blockSize), inBufferPt[1], blockSize * sizeof(float));
		// b - window the block in
		for(i = 0; i < sizeFFT; i++)
		{
			*(inFFTPt[1] + inputTimeR) = *(inShiftPt[1] + i) * *(analysisWindow + i);
			++inputTimeR;
			inputTimeR = inputTimeR & maskFFT;
		}
#if FFT_USED_FFTW
		rfftw_one(planIn, inFFTPt[1], inSpectralPt[1]);
#elif FFT_USED_APPLEVECLIB
	    vDSP_ctoz ((COMPLEX *)inFFTPt[1], 2, &split, 1, halfSizeFFT);
		vDSP_fft_zrip (setupReal, &split, 1, log2n, FFT_FORWARD);
		memcpy(inSpectralPt[1], split.realp, halfSizeFFT * sizeof(float));
		*(inSpectralPt[1] + halfSizeFFT) = *(split.imagp);
		for(i = halfSizeFFT + 1, j = halfSizeFFT - 1; i < sizeFFT; i++, j--)
			*(inSpectralPt[1] + i) = *(split.imagp + j);
#elif FFT_USED_INTEL
#endif
	}
	
	ProcessSpect();
	
	memcpy(outDisplaySpectPt[0], outSpectralPt[0], sizeFFT * sizeof(float));
	// d - IFFT
#if FFT_USED_FFTW
	rfftw_one(planOut, outSpectralPt[0], outShiftPt[0]);
#elif FFT_USED_APPLEVECLIB
	memcpy(split.realp, outSpectralPt[0], halfSizeFFT * sizeof(float));
	*(split.imagp) = *(outSpectralPt[0] + halfSizeFFT);
	for(i = halfSizeFFT + 1, j = halfSizeFFT - 1; i < sizeFFT; i++, j--)
		*(split.imagp + j) = *(outSpectralPt[0] + i);
	vDSP_fft_zrip (setupReal, &split, 1, log2n, FFT_INVERSE);
	vDSP_ztoc (&split, 1, ( COMPLEX * ) outShiftPt[0], 2, halfSizeFFT );
#elif FFT_USED_INTEL
#endif
	// e - overlap add
	for(i = 0; i < sizeFFT; i++)
	{
		outTemp = *(outShiftPt[0] + outputTimeL) * *(synthesisWindow + i);
		*(outBufferPt[0]+i) += outTemp;
		++outputTimeL;
		outputTimeL = outputTimeL & maskFFT;
	}
	if(channelMode == kStereoMode)
	{
		memcpy(outDisplaySpectPt[1], outSpectralPt[1], sizeFFT * sizeof(float));
		// d - IFFT
#if FFT_USED_FFTW
		rfftw_one(planOut, outSpectralPt[1], outShiftPt[1]);
#elif FFT_USED_APPLEVECLIB
		memcpy(split.realp, outSpectralPt[1], halfSizeFFT * sizeof(float));
		*(split.imagp) = *(outSpectralPt[1] + halfSizeFFT);
		for(i = halfSizeFFT + 1, j = halfSizeFFT - 1; i < sizeFFT; i++, j--)
			*(split.imagp + j) = *(outSpectralPt[1] + i);
		vDSP_fft_zrip (setupReal, &split, 1, log2n, FFT_INVERSE);
		vDSP_ztoc (&split, 1, ( COMPLEX * ) outShiftPt[1], 2, halfSizeFFT );
#elif FFT_USED_INTEL
#endif
		// e - overlap add
		for(i = 0; i < sizeFFT; i++)
		{
			outTemp = *(outShiftPt[1] + outputTimeR) * *(synthesisWindow + i);
			*(outBufferPt[1]+i) += outTemp;
			++outputTimeR;
			outputTimeR = outputTimeR & maskFFT;
		}
	}
}

void spectFlat::ProcessSpect(){
	// c - copy the spectra over - this is where we would typically do something
	if (channelMode == kStereoMode) {
		for (int i = 0; i < nMaxChannels; i++) {
			cartToPolar(inSpectralPt[i], polarSpectPt[i], decaySpectPt[i]);
			polarToCart(polarSpectPt[i], inSpectralPt[i], weightsPt[i], decaySpectPt[i]);
			outSpectralPt[i][0] = inSpectralPt[i][0];	// DC Component
			outSpectralPt[i][halfSizeFFT] = inSpectralPt[i][halfSizeFFT];	// Nyquist Frequency
			for (int j = 0; j < halfSizeFFT; j++) {
				outSpectralPt[i][j] = inSpectralPt[i][j];
				outSpectralPt[i][sizeFFT - j] = inSpectralPt[i][sizeFFT - j];
			}
		}
	}else if (channelMode == kMono2StereoMode){
		for (int i = 0; i < nMaxChannels; i++) {
			cartToPolar(inSpectralPt[0], polarSpectPt[0], decaySpectPt[0]);
			polarToCart(polarSpectPt[0], inSpectralPt[0], weightsPt[0], decaySpectPt[0]);
			outSpectralPt[i][0] = inSpectralPt[0][0];	// DC Component
			outSpectralPt[i][halfSizeFFT] = inSpectralPt[0][halfSizeFFT];	// Nyquist Frequency
			for (int j = 0; j < halfSizeFFT; j++) {
				outSpectralPt[i][j] = inSpectralPt[0][j];
				outSpectralPt[i][sizeFFT - j] = inSpectralPt[0][sizeFFT - j];
			}
		}
	}else{
		cartToPolar(inSpectralPt[0], polarSpectPt[0], decaySpectPt[0]);
		polarToCart(polarSpectPt[0], inSpectralPt[0], weightsPt[0], decaySpectPt[0]);
		outSpectralPt[0][0] = inSpectralPt[0][0];	// DC Component
		outSpectralPt[0][halfSizeFFT] = inSpectralPt[0][halfSizeFFT];	// Nyquist Frequency
		for (int j = 0; j < halfSizeFFT; j++) {
			outSpectralPt[0][j] = inSpectralPt[0][j];
			outSpectralPt[0][sizeFFT - j] = inSpectralPt[0][sizeFFT - j];
		}
	}
}

void spectFlat::initHammingWindows(void){
	long	index;
	float	a = 0.54f, b	= 0.46f;
	winType = kHamming;
	// a - make two hamming windows
	for (index = 0; index < sizeFFT; index++)
		synthesisWindow[index] = analysisWindow[index] = a - b*cosf(twoPi*index/(sizeFFT - 1));
}

void spectFlat::scaleWindows(void){
	long	index;
	float	sum, analFactor, synthFactor;
	

	// b - scale the windows
	sum = 0.0f;
	for (index = 0; index < sizeFFT; index++)
		sum += analysisWindow[index];
	
	synthFactor = analFactor = 2.0f/sum;
	
    for (index = 0; index < sizeFFT; index++)
    {
		analysisWindow[index] *= analFactor;
		synthesisWindow[index] *= synthFactor;
	}
	sum = 0.0;
	for (index = 0; index < sizeFFT; index += blockSize)
   		sum += synthesisWindow[index]*synthesisWindow[index];
   	// i think this scaling factor is only needed for vector lib
#if FFT_USED_FFTW
   	sum = 1.0f/(sum*sizeFFT);
#elif FFT_USED_APPLEVECLIB
   	sum = 1.0f/(sum*sizeFFT*2.0f);
#elif FFT_USED_INTEL
   	sum = 1.0f/(sum*sizeFFT);
#endif
	for (index = 0; index < sizeFFT; index++)
   		synthesisWindow[index] *= sum;
}

void spectFlat::cartToPolar(float *spectrum, float *polarSpectrum, float *decay){
	
	long	realIndex, imagIndex, ampIndex, phaseIndex;
	float	realPart, imagPart;
	long	bandNumber;
	magSum = 0.0f;
	
	for (bandNumber = 0; bandNumber <= halfSizeFFT; bandNumber++)
	{
		realIndex = bandNumber;
		imagIndex = sizeFFT - bandNumber;
		ampIndex = bandNumber<<1;
		phaseIndex = ampIndex + 1;
		if(bandNumber == 0)
		{
			realPart = spectrum[0];
			imagPart = 0.0;
		}
		else if(bandNumber == halfSizeFFT)
		{
			realPart = spectrum[halfSizeFFT];
			imagPart = 0.0;
		}
		else
		{
			realPart = spectrum[realIndex];
			imagPart = spectrum[imagIndex];
		}
		/*
		 * compute magnitude & phase value from real and imaginary parts
		 */
		
		polarSpectrum[ampIndex] = hypot(realPart, imagPart);
		//
		//
		if(polarSpectrum[ampIndex] < 0.0000001f)
			polarSpectrum[phaseIndex] = 0.0;
		else
			polarSpectrum[phaseIndex] = atan2f(imagPart, realPart);
		float f = polarSpectrum[ampIndex];
		decay[bandNumber] = decay[bandNumber] + (f - decay[bandNumber]/decaySize);
		//magSum += polarSpectrum[ampIndex];
		magSum += decay[bandNumber];
	}
	magSum *= oneOverFFTSize;
}


void spectFlat::polarToCart(float *polarSpectrum, float *spectrum, float* weight, float *decay){
	
	float	realValue, imagValue;
	long	bandNumber, realIndex, imagIndex, ampIndex, phaseIndex;
	
	/*
	 * convert to cartesian coordinates, complex pairs
	 */
    for (bandNumber = 0; bandNumber <= halfSizeFFT; bandNumber++)
    {
		realIndex = bandNumber;
		imagIndex = sizeFFT - bandNumber;
		ampIndex = bandNumber<<1;
		phaseIndex = ampIndex + 1;
				
		float tmpW;
		if (magSum == 0 || decay[bandNumber] == 0) {
			tmpW = 1.0f;
		}else{
			tmpW = magSum/decay[bandNumber];
		}
		
		weight[bandNumber] += stepSize*(tmpW - weight[bandNumber]);
		
		// Not crashing Ableton but the "not-white" == pink noise part is not working as expected.
		if (white) {
			polarSpectrum[ampIndex] = polarSpectrum[ampIndex] + max*(weight[bandNumber]-1.0f)*polarSpectrum[ampIndex];
		}else{
			polarSpectrum[ampIndex] = polarSpectrum[ampIndex] + max*(1.0f/((float)bandNumber+1.0f))*(weight[bandNumber]-1.0f)*polarSpectrum[ampIndex];
		}
		// -------------------------------------------------------------------------
		
		if (polarSpectrum[ampIndex]<0.0000001f) {
			polarSpectrum[ampIndex] = 0;
		}

		if(polarSpectrum[ampIndex] == 0.0)
		{
			realValue = 0.0;
			imagValue = 0.0;
		}
		else if(bandNumber == 0 || bandNumber == halfSizeFFT)
    	{
			realValue = polarSpectrum[ampIndex] * cosf(polarSpectrum[phaseIndex]);
    		imagValue = 0.0;
    	}
    	else
		{
			realValue = polarSpectrum[ampIndex] * cosf(polarSpectrum[phaseIndex]);
			imagValue = polarSpectrum[ampIndex] * sinf(polarSpectrum[phaseIndex]);
		}
		
		if(bandNumber == halfSizeFFT)
			realIndex = halfSizeFFT;
		spectrum[realIndex] = realValue;
		if(bandNumber != halfSizeFFT && bandNumber != 0)
			spectrum[imagIndex] = imagValue;
	}
	memset(polarSpectrum, 0, kSizeFFT*sizeof(float));
}



