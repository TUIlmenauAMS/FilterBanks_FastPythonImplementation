# Fast Filter Bank Implementations using Python 

Programs for the implementation of fast filter banks, Low Delay filter banks, Modified Discrete Cosine Transform filter banks, and Pseudo Quadrature Mirror filter banks.

The advantages of the fast imnplementation are: 
* fastest implementation without the need of specialized hardware
* low memory requirement for filtering and signal processing
* real time processing, one signal block goes in, one signal block comes out

For the theoretical background see the book:

Gerald Schuller: "Filter Banks and Audio Coding - Compressing Audio Signals Using Python"
Springer 2020, 
ISBN: 978-3-030-51249-1 (e-book)
ISBN: 978-3-030-51251-4 (softcover or hardcover)

The programs run under both, Python2 and Python3.
Among others they need OpenCV (cv2) for the live waterfall spectrogram.

Under Python2 it is installed with:

sudo apt install python-opencv

in Python3 it is installed with:

sudo apt install python3-pip

sudo pip3 install opencv-python


