# Fast Pseudo Quadrature Mirror Filter Bank Filter Bank Implementation 

Programs for the implementation of PQMF  filter banks.
For a background on PQMF  filter banks, see e.g.:
https://ccrma.stanford.edu/~jos/sasp/Pseudo_QMF_Cosine_Modulation_Filter.html
the Book:  Spanias, Painter: “Audio Signal Processing and Coding”, Wiley Press
and our lecture [Audio Coding] 
(https://www.tu-ilmenau.de/mt/lehrveranstaltungen/lehre-fuer-master-mt/audio-coding/)
slides lecture 07.

PQMF filter banks have longer filters than MDCT filter banks with the same number of subbands, hence have (and need) a higher stopbband attenuation, but at the expense of a higher system delay (over the analysis and synthesis filter bank), and have only a near perfect reconstruction, which means there is always a reconstruction error, although hopefully very small (for that the high stopband attenuation is needed).

The advantages of the fast implementation are: 
* fastest implementation without the need of specialized hardware
* low memory requirement for filtering and signal processing
* real time processing, one signal block goes in, one signal block comes out


## Getting Started
The program "qmf_realtime.py" contains the fast computation of the PQMF analysis and synthesis filter bank.
In the beginning, it sets the number of subbands N=64, as is usual for instance in MPEG coders, and in the "main" section it reads in the prototype filter coefficients from file 'qmf.dat'.

Execute it with:
python qmf_realtime.py

First it runs a test of the filter bank by inputting a single pulse into the synthesis filter bank filter 0, plots the resulting impulse response, and its corresponding frequency response.
After closing these windows, it takes the sound file "test.wav" as input, then it runs the analysis and synthesis filter bank on its signal in a for loop, and writes the reconstructed signal into the file "testrek.wav", until the entire signal is processed.


Gerald Schuller, Februrary 2018.

