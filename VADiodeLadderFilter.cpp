/*
	RackAFX(TM)
	Applications Programming Interface
	Derived Class Object Implementation
*/


#include "VADiodeLadderFilter.h"


/* constructor()
	You can initialize variables here.
	You can also allocate memory here as long is it does not
	require the plugin to be fully instantiated. If so, allocate in init()

*/
CVADiodeLadderFilter::CVADiodeLadderFilter()
{
	
	m_dGAMMA = 0.0;
	m_dK = 0.0;

	// our feedback S values (global)
	m_dSG1 = 0.0; 
	m_dSG2 = 0.0;
	m_dSG3 = 0.0; 
	m_dSG4 = 0.0; 

	// Filter coeffs that are constant
	// set a0s
	m_LPF1.m_da0 = 1.0;
	m_LPF2.m_da0 = 0.5;
	m_LPF3.m_da0 = 0.5;
	m_LPF4.m_da0 = 0.5;

	// last LPF has no feedback path
	m_LPF4.m_dGamma = 1.0;
	m_LPF4.m_dDelta = 0.0;
	m_LPF4.m_dEpsilon = 0.0;
	m_LPF4.setFeedback(0.0);

}


/* destructor()
	Destroy variables allocated in the contructor()

*/
CVADiodeLadderFilter::~CVADiodeLadderFilter(void)
{


}

/*
initialize()
	Called by the client after creation; the parent window handle is now valid
	so you can use the Plug-In -> Host functions here (eg sendUpdateUI())
	See the website www.willpirkle.com for more details
*/
bool CVADiodeLadderFilter::initialize()
{
	// Add your code here

	return true;
}



/* prepareForPlay()
	Called by the client after Play() is initiated but before audio streams

	You can perform buffer flushes and per-run intializations.
	You can check the following variables and use them if needed:

	m_nNumWAVEChannels;
	m_nSampleRate;
	m_nBitDepth;

	NOTE: the above values are only valid during prepareForPlay() and
		  processAudioFrame() because the user might change to another wave file,
		  or use the sound card, oscillators, or impulse response mechanisms

    NOTE: if you allocte memory in this function, destroy it in ::destroy() above
*/
bool __stdcall CVADiodeLadderFilter::prepareForPlay()
{
	// Add your code here:

	// this flushes all storage registers in filters including feedback register
	reset();

	// set the initial coeffs
	updateFilter();

	return true;
}

void CVADiodeLadderFilter::updateFilter()
{
	// calculate alphas
	double wd = 2*pi*m_dFc;          
	double T  = 1/(float)m_nSampleRate;             
	double wa = (2/T)*tan(wd*T/2); 
	double g = wa*T/2;  

	// Big G's
	double G1, G2, G3, G4;

	G4 = 0.5*g/(1.0 + g);
	G3 = 0.5*g/(1.0 + g - 0.5*g*G4);
	G2 = 0.5*g/(1.0 + g - 0.5*g*G3);
	G1 = g/(1.0 + g - g*G2);
	
	// our big G value GAMMA
	m_dGAMMA = G4*G3*G2*G1;
	
	m_dSG1 =  G4*G3*G2; 
	m_dSG2 =  G4*G3; 
	m_dSG3 =  G4; 
	m_dSG4 =  1.0; 

	// set alphas
	m_LPF1.m_dAlpha = g/(1.0 + g);
	m_LPF2.m_dAlpha = g/(1.0 + g);
	m_LPF3.m_dAlpha = g/(1.0 + g);
	m_LPF4.m_dAlpha = g/(1.0 + g);

	// set betas
	m_LPF1.m_dBeta = 1.0/(1.0 + g - g*G2);
	m_LPF2.m_dBeta = 1.0/(1.0 + g - 0.5*g*G3);
	m_LPF3.m_dBeta = 1.0/(1.0 + g - 0.5*g*G4);
	m_LPF4.m_dBeta = 1.0/(1.0 + g);
	
	// set gammas
	m_LPF1.m_dGamma = 1.0 + G1*G2;
	m_LPF2.m_dGamma = 1.0 + G2*G3;
	m_LPF3.m_dGamma = 1.0 + G3*G4;
	// m_LPF4.m_dGamma = 1.0; // constant - done in constructor
	
	// set deltas
	m_LPF1.m_dDelta = g;
	m_LPF2.m_dDelta = 0.5*g;
	m_LPF3.m_dDelta = 0.5*g;
	// m_LPF4.m_dDelta = 0.0; // constant - done in constructor

	// set epsilons
	m_LPF1.m_dEpsilon = G2;
	m_LPF2.m_dEpsilon = G3;
	m_LPF3.m_dEpsilon = G4;
	// m_LPF4.m_dEpsilon = 0.0; // constant - done in constructor
}

// do the filter
double CVADiodeLadderFilter::doFilter(double xn)
{
//	m_LPF4.setFeedback(0.0); // constant - done in constructor
	m_LPF3.setFeedback(m_LPF4.getFeedbackOutput());
	m_LPF2.setFeedback(m_LPF3.getFeedbackOutput());
	m_LPF1.setFeedback(m_LPF2.getFeedbackOutput());

	// form input
	double SIGMA = m_dSG1*m_LPF1.getFeedbackOutput() + 
				   m_dSG2*m_LPF2.getFeedbackOutput() +
				   m_dSG3*m_LPF3.getFeedbackOutput() +
				   m_dSG4*m_LPF4.getFeedbackOutput();

	// "cheap" nonlinear model; just process input
	if(m_NonLinearProcessing == ON)
	{
		// the 1/tanh(sat) is to normalize the function so x = [-1..+1] --> y = [-1..+1]
		// Normalized Version
		if(m_uNLPType == NORM)
			xn = (1.0/tanh(m_dSaturation))*tanh(m_dSaturation*xn);
		else
			xn = tanh(m_dSaturation*xn);
	}

	// form the input to the loop
	double un = (xn - m_dK*SIGMA)/(1 + m_dK*m_dGAMMA);

	// cascade of series filters
	return m_LPF4.doFilter(m_LPF3.doFilter(m_LPF2.doFilter(m_LPF1.doFilter(un))));
}


/* processAudioFrame

// ALL VALUES IN AND OUT ON THE RANGE OF -1.0 TO + 1.0

LEFT INPUT = pInputBuffer[0];
RIGHT INPUT = pInputBuffer[1]

LEFT INPUT = pInputBuffer[0]
LEFT OUTPUT = pOutputBuffer[1]

*/
bool __stdcall CVADiodeLadderFilter::processAudioFrame(float* pInputBuffer, float* pOutputBuffer, UINT uNumInputChannels, UINT uNumOutputChannels)
{
	// Monophonic filter
	//
	// Do LEFT (MONO) Channel; there is always at least one input/one output
	// (INSERT Effect)
	pOutputBuffer[0] = doFilter((double)pInputBuffer[0]);

	// Mono-In, Stereo-Out (AUX Effect)
	if(uNumInputChannels == 1 && uNumOutputChannels == 2)
		pOutputBuffer[1] = pOutputBuffer[0];

	// Stereo-In, Stereo-Out (INSERT Effect)
	if(uNumInputChannels == 2 && uNumOutputChannels == 2)
		pOutputBuffer[1] = pOutputBuffer[0];


	return true;
}



//


/* processAudioBuffer

	// ALL VALUES IN AND OUT ON THE RANGE OF -1.0 TO + 1.0

	The I/O buffers are interleaved depending on the number of channels. If uNumChannels = 2, then the
	buffer is L/R/L/R/L/R etc...

	if uNumChannels = 6 then the buffer is L/R/C/Sub/BL/BR etc...

	It is up to you to decode and de-interleave the data.

	To use this function set m_bWantBuffers = true in your constructor.

	******************************
	********* IMPORTANT! *********
	******************************
	If you are going to ultimately make this a VST Compatible Plug-In and you want to process
	buffers, you need to override the NEXT function below:

	processVSTAudioBuffer()


	This function (processRackAFXAudioBuffer) is not supported in the VST wrapper because
	the VST buffer sizes no maximum value. This would require the use of dynamic buffering
	in the callback which is not acceptable for performance!
*/
bool  CVADiodeLadderFilter::processRackAFXAudioBuffer(float* pInputBuffer, float* pOutputBuffer,
													   UINT uNumInputChannels, UINT uNumOutputChannels,
													   UINT uBufferSize)
{

	for(UINT i=0; i<uBufferSize; i++)
	{
		// pass through code
		pOutputBuffer[i] = pInputBuffer[i];
	}


	return true;
}



/* processVSTAudioBuffer

	// ALL VALUES IN AND OUT ON THE RANGE OF -1.0 TO + 1.0

	NOTE: You do not have to implement this function if you don't want to; the processAudioFrame()
	will still work; however this using function will be more CPU efficient for your plug-in, and will
	override processAudioFrame().

	To use this function set m_bWantVSTBuffers = true in your constructor.

	The VST input and output buffers are pointers-to-pointers. The pp buffers are the same depth as uNumChannels, so
	if uNumChannels = 2, then ppInputs would contain two pointers,

		inBuffer[0] = a pointer to the LEFT buffer of data
		inBuffer[1] = a pointer to the RIGHT buffer of data

	Similarly, outBuffer would have 2 pointers, one for left and one for right.

	For 5.1 audio you would get 6 pointers in each buffer.

*/
//bool __stdcall CVADiodeLadderFilter::processVSTAudioBuffer(float** inBuffer, float** outBuffer, UINT uNumChannels, int inFramesToProcess)
//{
//	// PASS Through example
//	// MONO First
//	float* pInputL  = inBuffer[0];
//	float* pOutputL = outBuffer[0];
//	float* pInputR  = NULL;
//	float* pOutputR = NULL;
//
//	// if STEREO,
//	if(inBuffer[1])
//		pInputR = inBuffer[1];
//
//	if(outBuffer[1])
//		pOutputR = outBuffer[1];
//
//	// Process audio by de-referencing ptrs
//	// this is siple pass through code
//	while (--inFramesToProcess >= 0)
//	{
//		// Left channel processing
//		*pOutputL = *pInputL;
//
//		// If there is a right channel
//		if(pInputR)
//			*pOutputR = *pInputR;
//
//		// advance pointers
//		pInputL++;
//		pOutputL++;
//		if(pInputR) pInputR++;
//		if(pOutputR) pOutputR++;
//	}
//	// all OK
//	return true;
//}



























