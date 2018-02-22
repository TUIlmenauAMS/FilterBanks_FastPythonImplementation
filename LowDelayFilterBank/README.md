# Fast Low Delay Filter Bank Implementation 

Programs for the implementation of Low Delay filter banks.
For background on Low Delay filter banks see:
https://www.idmt.fraunhofer.de/content/dam/idmt/documents/IL/Personal%20Websites/Schuller/publications/tsp8-96.pdf
and
https://www.idmt.fraunhofer.de/content/dam/idmt/documents/IL/Personal%20Websites/Schuller/publications/tsp3-00.pdf
and our lecture [Multirate Signal Processing] 
(https://www.tu-ilmenau.de/mt/lehrveranstaltungen/lehre-fuer-master-mt/multirate-signal-processing/)
slides lecture 15.

The advantages of the fast implementation are: 
* fastest implementation without the need of specialized hardware
* low memory requirement for filtering and signal processing
* real time processing, one signal block goes in, one signal block comes out


## Getting Started
The program "pyrecplayfastLDFB.py" contains the fast computation of the LDFB analysis and synthesis filter bank.
In the "main" section, it sets the number of subbands N=512.

Execute it with:
python pyrecplayfastLDFB.py

First it runs a test of the filter bank by inputting a single pulse into the synthesis filter banks filter 0, plots the resulting impulse response, and its resulting frequency response.
After closing these windows, it starts a window with a live waterfall spectrogram. It shows the 512 subbands of the LDFB analysis filter bank from the sound card microphone input. 
Then the subband signals are fed into the synthesis filter bank to reconstruct the signal, and fed into the sound card speaker. Reduce the volume if necessary to avoid a feedback loop.
The live spectrogram is closed by pressing the key "q" while the spectrogram window is active (selected by the mouse).

The program "fastLDFBstereo_anasyn.py" takes a stereo sound file as input, then it runs the analysis and synthesis filter bank on its signal in a for loop, until the entire signal is processed. Then it play back the reconstructed signal, and displays the time needed for processing it.

To start, simply execute:
python fastLDFBstereo_anasyn.py


Gerald Schuller, Februrary 2018.

