# specFlat
An AudioUnit inspired by the concept of spectral flatness. 

# author & contact
Cheng-i Wang, chw160@ucsd.edu

# description
SpecFlat is an Audio Unit developed with JUCE. The plug-in transforms the spectrum of the audio input to maximize/minimize the spectral flatness measure. Maximizing the spectral flatness is the same as whitening the spectrum, while minimize is not well-defined given an arbitrary input, the minimization tends to reinforce the stronger spectral part in the spectrum. 
