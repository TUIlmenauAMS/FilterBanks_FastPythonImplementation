# Fast MDCT Filter Bank Implementation using Python 

Programs for the implementation of MDCT filter banks.
For background on MDCT filter banks, see e.g. our lecture [Multirate Signal Processing] 
(https://www.tu-ilmenau.de/mt/lehrveranstaltungen/lehre-fuer-master-mt/multirate-signal-processing/)
slides lecture 14.

The advantages of the fast imnplementation are: 
* fastest implementation without the need of specialized hardware
* low memory requirement for filtering and signal processing
* real time processing, one signal block goes in, one signal block comes out


## Getting Started
The program "pyrecplayfastMDCT.py" contains the fast computation of the MDCT analysis and synthesis filter bank.
In the "main" section, it sets the number of subbands N=512.
It can be set to any even number, and this produces a "sine window" of length 2N.
Execute it with:
python pyrecplayfastMDCT.py

First it runs a test of the filter bank by inputting a single pulse into the synthesis filter bank, filter 0, plots the resulting impulse response, and its corresponding frequency response.
After closing these windows, it starts a window with a live waterfall spectrogram. It shows the 512 subbands of the MDCT analysis filter bank from the sound card microphone input (each column is a subband). 
Then the subband signals are fed into the synthesis filter bank to reconstruct the signal, and fed into the sound card speaker. Reduce the volume if necessary to avoid a feedback loop.
The live spectrogram is closed by pressing the key "q" while the spectrogram window is active (selected by the mouse).

The program "fastMDCTanasyn.py" takes a sound file as input, then it runs the analysis and synthesis filter bank on its signal in a for loop, until the entire signal is processed. Then it plays back the reconstructed signal, and displays the time needed for processing it.

To start, simply execute:
python fastMDCTanasyn.py


Gerald Schuller, Februrary 2018.

