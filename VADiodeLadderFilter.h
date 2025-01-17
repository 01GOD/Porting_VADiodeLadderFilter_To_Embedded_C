

#pragma once

// base class
//#include "plugin.h"
#include "VAOnePoleFilterEx.h"

// abstract base class for RackAFX filters
class CVADiodeLadderFilter
{
public:
//	// RackAFX Plug-In API Member Methods:
//	// The followung 5 methods must be impelemented for a meaningful Plug-In
//	//
//	// 1. One Time Initialization
//	CVADiodeLadderFilter();
//
//	// 2. One Time Destruction
//	virtual ~CVADiodeLadderFilter(void);
//
//	// 3. The Prepare For Play Function is called just before audio streams
//	virtual bool __stdcall prepareForPlay();

	// 4. processAudioFrame() processes an audio input to create an audio output
	virtual bool processAudioFrame(float* pInputBuffer, float* pOutputBuffer, UINT uNumInputChannels, UINT uNumOutputChannels);

//	// 5. userInterfaceChange() occurs when the user moves a control.
//	virtual bool __stdcall userInterfaceChange(int nControlIndex);
//
//
	// OPTIONAL ADVANCED METHODS ------------------------------------------------------------------------------------------------
	// These are more advanced; see the website for more details
	//
	// 6. initialize() is called once just after creation; if you need to use Plug-In -> Host methods
	//				   such as sendUpdateGUI(), you must do them here and NOT in the constructor
//	virtual bool __stdcall initialize();
//
//	// 7. joystickControlChange() occurs when the user moves a control.
//	virtual bool __stdcall joystickControlChange(float fControlA, float fControlB, float fControlC, float fControlD, float fACMix, float fBDMix);
//
//	// 8. process buffers instead of Frames:
//	// NOTE: set m_bWantBuffers = true to use this function
//	virtual bool __stdcall processRackAFXAudioBuffer(float* pInputBuffer, float* pOutputBuffer, UINT uNumInputChannels, UINT uNumOutputChannels, UINT uBufferSize);
//
//	// 9. rocess buffers instead of Frames:
//	// NOTE: set m_bWantVSTBuffers = true to use this function
//	virtual bool __stdcall processVSTAudioBuffer(float** inBuffer, float** outBuffer, UINT uNumChannels, int inFramesToProcess);
//
//	// 10. MIDI Note On Event
//	virtual bool __stdcall midiNoteOn(UINT uChannel, UINT uMIDINote, UINT uVelocity);
//
//	// 11. MIDI Note Off Event
//	virtual bool __stdcall midiNoteOff(UINT uChannel, UINT uMIDINote, UINT uVelocity, bool bAllNotesOff);
//
//
//	// 12. MIDI Modulation Wheel uModValue = 0 -> 127
//	virtual bool __stdcall midiModWheel(UINT uChannel, UINT uModValue);
//
//	// 13. MIDI Pitch Bend
//	//					nActualPitchBendValue = -8192 -> 8191, 0 is center, corresponding to the 14-bit MIDI value
//	//					fNormalizedPitchBendValue = -1.0 -> +1.0, 0 is at center by using only -8191 -> +8191
//	virtual bool __stdcall midiPitchBend(UINT uChannel, int nActualPitchBendValue, float fNormalizedPitchBendValue);
//
//	// 14. MIDI Timing Clock (Sunk to BPM) function called once per clock
//	virtual bool __stdcall midiClock();
//
//
//	// 15. all MIDI messages -
//	// NOTE: set m_bWantAllMIDIMessages true to get everything else (other than note on/off)
//	virtual bool __stdcall midiMessage(unsigned char cChannel, unsigned char cStatus, unsigned char cData1, unsigned char cData2);
//
//	// 16. initUI() is called only once from the constructor; you do not need to write or call it. Do NOT modify this function
//	virtual bool __stdcall initUI();
//


	// Add your code here: ----------------------------------------------------------- //
	CVAOnePoleFilterEx m_LPF1;
	CVAOnePoleFilterEx m_LPF2;
	CVAOnePoleFilterEx m_LPF3;
	CVAOnePoleFilterEx m_LPF4;

	double m_dGAMMA; // Gamma see App Note

	// our feedback S values (global)
	double m_dSG1; 
	double m_dSG2; 
	double m_dSG3; 
	double m_dSG4; 

	void reset()
	{
		m_LPF1.reset(); m_LPF2.reset();
		m_LPF3.reset(); m_LPF4.reset();

		m_LPF1.setFeedback(0.0); m_LPF2.setFeedback(0.0); 
		m_LPF3.setFeedback(0.0); m_LPF4.setFeedback(0.0); 
	}

	// recalc the coeffs
	void updateFilter();
	
	// do the filter
	double doFilter(double xn);
	// END OF USER CODE -------------------------------------------------------------- //


	// ADDED BY RACKAFX -- DO NOT EDIT THIS CODE!!! ----------------------------------- //
	//  **--0x07FD--**

	double m_dFc;
	double m_dK;
	double m_dSaturation;
	UINT m_NonLinearProcessing;
	enum{OFF,ON};
	UINT m_uNLPType;
	enum{NORM,REG};

	// **--0x1A7F--**
	// ------------------------------------------------------------------------------- //

};














