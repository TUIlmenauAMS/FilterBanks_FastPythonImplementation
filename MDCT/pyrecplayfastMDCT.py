# coding: utf-8 
#Program for the fast implementation of the MDCT and Low Delay analysis and synthesis filter bank.
#Gerald Schuller, Nov. 2014.
#Algorithm according to:
#G. Schuller and T. Karp: "Modulated Filter Banks with Arbitrary System Delay: Efficient Implementations and #the Time-Varying Case", IEEE Transactions on Signal Processing, March 2000, pp. 737–748 
#updated July 2016

import pyaudio
import struct
#import math
#import array
import numpy as np
#import sys
#import wave
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
#import pylab
import cv2
import scipy.fftpack as spfft
from scipy.signal import freqz
import time


#The Low Delay Filter Bank:------------------------------



def initFB(N):
   #initialize memory for the number of subbands N:
   Dmatrix.z=np.zeros(int(N/2))
   Dinvmatrix.z=np.zeros(int(N/2))
   Gmatrix.z=np.zeros((2,int(N/2)))
   Ginvmatrix.z=np.zeros((2,int(N/2)))
   H2matrix.z=np.zeros(N)
   H2invmatrix.z=np.zeros(N)
   return

#The D(z) matrix:
def Dmatrix(samples):
#implementation of the delay matrix D(z)
   #Delay elements:
   N=len(samples)
   out=np.zeros(N)
   out[0:int(N/2)]=Dmatrix.z
   Dmatrix.z=samples[0:int(N/2)]
   out[int(N/2):N]=samples[int(N/2):N]
   return out




#The inverse D(z) matrix:
def Dinvmatrix(samples):
#implementation of the delay matrix D(z)
   #Delay elements:
   N=len(samples)
   out=np.zeros(N)
   out[int(N/2):N]=Dinvmatrix.z
   Dinvmatrix.z=samples[int(N/2):N]
   out[0:int(N/2)]=samples[0:int(N/2)]
   return out





#The symmetric F Matrix function:
#implements multiplication samples*symFmatrix 
def symFmatrix(samples,fcoeff):
   sym=1.0;
   N=len(samples)
   out=np.zeros(N)
   out[0:int(N/2)]=(fcoeff[0:int(N/2)]*samples[0:int(N/2)])[::-1] +fcoeff[int(N/2):N]*samples[int(N/2):N]
   out[int(N/2):N]=(fcoeff[N:int(N+N/2)]*samples[0:int(N/2)]) 
   #+-1=det([a,b;c,d]) =a*d-b*c => d=(+-1+b*c)/a:
   ff= (-sym*np.ones(int(N/2))+ fcoeff[N:int(1.5*N)][::-1]*fcoeff[int(N/2):N])/fcoeff[0:int(N/2)][::-1]
   out[int(N/2):N]= out[int(N/2):N] +(ff*samples[int(N/2):N])[::-1]
   return out

def symFinvmatrix(samples,fcoeff):
   #inverse symFamtrix, uses order for the synthesis F matrix as shown in Audio Coding lecture FB2,
   #but with coefficients in reverse order, since the lecture has h as a window and not impulse response.
   #That is also why the negative signa has to be moved from the end to the beginning:
   sym=1.0
   N=len(samples)
   out=np.zeros(N)
   ff= (-sym*np.ones(int(N/2))+ fcoeff[N:int(1.5*N)][::-1]*fcoeff[int(N/2):N])/fcoeff[0:int(N/2)][::-1]

   out[0:int(N/2)]=-ff[::-1]*(samples[0:int(N/2)][::-1])+fcoeff[int(0.5*N):N][::-1]*samples[int(N/2):N]
   out[int(N/2):N]=fcoeff[N:int(1.5*N)][::-1]*samples[0:int(N/2)]- fcoeff[0:int(0.5*N)][::-1]*(samples[int(N/2):N][::-1])
   return out

#fcoeff=np.sin(np.pi/(2*N)*(np.arange(0,2*N)+0.5))
#Fmat=np.zeros((N,N))
#Fmat[0:N/2,0:N/2]=np.fliplr(np.diag(fcoeff[0:N/2]))
#Fmat[N/2:N,0:N/2]=np.diag(fcoeff[N/2:N])
#Fmat[0:N/2,N/2:N]=np.diag(fcoeff[N:(N+N/2)])
#Fmat[N/2:N,N/2:N]=-np.fliplr(np.diag(fcoeff[(N+N/2):(2*N)]))

#The inverse F matrix:
#Finv=np.linalg.inv(Fmat)


#The G_i(z) matrix:
def Gmatrix(i,samples,ecoeff):
#implementation of the delay matrix G(z)
   N=len(samples)
   #Anti-diagonal ones, flip input:
   out=samples[::-1]
   #Delay elements and coeff in the upper half of the diagonal:
   
   out[0:int(N/2)]=out[0:int(N/2)]+Gmatrix.z[i] * ecoeff
   Gmatrix.z[i]=samples[0:int(N/2)]
   return out



#The inverse G_i(z) matrix:
def Ginvmatrix(i,samples,ecoeff):
#implementation of the delay matrix G(z)
   N=len(samples)
   #Anti-diagonal ones, flip input:
   out=samples[::-1]
   #Delay elements and flipped, neg. coeff in the lower half of the diagonal:
   out[int(N/2):N]=out[int(N/2):N] - Ginvmatrix.z[i] * ecoeff[::-1]
   Ginvmatrix.z[i]=samples[int(N/2):N]
   return out



#The H2(z) matrix:
def H2matrix(samples,H2coeff):
#implementation of the delay matrix G(z)
   N=len(samples)
   #Anti-diagonal delays, flip delayed input:
   out=H2matrix.z[::-1]
   #input mult. with coeff in the upper half of the diagonal:
   out[0:int(N/2)]=out[0:int(N/2)]+ samples[0:int(N/2)]* H2coeff
   H2matrix.z=samples
   return out



#The inverse H2(z) matrix:
def H2invmatrix(samples,H2coeff):
#implementation of the delay matrix G(z)
   N=len(samples)
   #Anti-diagonal delays, flip delayed input:
   out=H2invmatrix.z[::-1]
   #input mult. with neg. flipped coeff in the lower half of the diagonal:
   
   out[int(N/2):N]=out[int(N/2):N]- samples[int(N/2):N]* H2coeff[::-1]
   H2invmatrix.z=samples
   return out




#The DCT4 transform:
def DCT4(samples):
   #use a DCT3 to implement a DCT4:
   N=len(samples)
   samplesup=np.zeros(2*N)
   #upsample signal:
   samplesup[1::2]=samples
   y=spfft.dct(samplesup,type=3)/2
   return y[0:N]


#The complete LDFB, Analysis:
def LDFB(samples,fb):
   #samples: N samples of input signal
   #fb: filter bank coefficients
   N=len(samples)
   #load LDFB coefficients:
   #Fmatrix:
   fcoeff=fb[0:int(1.5*N)]
   #G_0 matrix:
   #g0coeff=fb[(1.5*N):(2*N)]
   #G1 matrix:
   #g1coeff=fb[(2*N):(2.5*N)]
 
   y=symFmatrix(samples, fcoeff)
   y=Dmatrix(y)
   #y=Gmatrix(0,y,g0coeff)
   #y=Gmatrix(1,y,g1coeff)
   y=DCT4(y)
   return y

#The inverse LDFB, synthesis:
def LDFBinv(y,fb): 
   #y: N subband samples or frequency coefficients
   #fb: filter bank coefficients
   N=len(y)
   #load LDFB coefficients:
   #Fmatrix:
   fcoeff=fb[0:int(1.5*N)]
   #G_0 matrix:
   #g0coeff=fb[(1.5*N):(2*N)]
   #G1 matrix:
   #g1coeff=fb[(2*N):(2.5*N)]
    
   #inverse DCT4 is identical to DCT4:
   x=DCT4(y)*2/N
   #inverse D(z) matrix:
   x=Dinvmatrix(x)
   #inverse F matrix:
   x=symFinvmatrix(x,fcoeff)
   return x

if __name__ == '__main__':
        #---------------------------
        #Test Filter Bank:
        """
        print("test symFmatrix:")
        x=np.arange(1,9)
        N=8;
        initFB(N)
        fcoeff=np.arange(1,13)
        xF=symFmatrix(x,fcoeff)
        #print(xF)
        #x=np.zeros(8)
        #x[7]=1.0
        xrek=symFinvmatrix(xF,fcoeff)
        print(xrek)

        N=1024
        initFB(N)
        fb=np.sin(np.pi/(2*N)*(np.arange(0,1.5*N)+0.5))
        print("Test FB")

        #Test Perfect reconstruction:
        x=np.arange(10*N)
        xrek=np.zeros(10*N)
        for m in range(10):
          print(m)
          a=time.time()
          y=LDFB(x[(m*N):((m+1)*N)],fb)
          b=time.time()
          xrek[(m*N):((m+1)*N)]=LDFBinv(y,fb)
          c=time.time()
          print("Analsysis time:", b-a)
          print("Synthesis time:", c-b)
        #Plots reconstructed ramp if ok:
        fig1 = plt.figure()
        plt.plot(xrek)
        fig1.show()

        #Test impulse response in synthesis:
        #Set memories of synthesis to zero:
        H2invmatrix.z=np.zeros(N)
        Ginvmatrix.z[0]=np.zeros(N/2)
        Ginvmatrix.z[1]=np.zeros(N/2)
        Dinvmatrix.z=np.zeros(N/2)
        #Set impulse in subband 0 for synthesis:
        y=np.ones((N,10))*0.0;
        y[0,0]=1.0
        xrek=np.ones(10*N)*0.0
        for m in range(10):
           xrek[(m*N):((m+1)*N)]=LDFBinv(y[:,m])
        #Plots synthesis impulse response of subband 0:
        fig2 = plt.figure()
        plt.plot(xrek)
        fig2.show()
        """
        #Runs real time audio processing:----------------------------

        #Load the LDFB coefficients from text file:
        #256 taps, 127 delay, 64 subbands:
        #fb=np.loadtxt('fbsy256t127d64bbitb.mat')
        

        N=512 #Number of subbands and block size
        display=True
        #example: Sine window:
        fb=np.sin(np.pi/(2*N)*(np.arange(0,int(1.5*N))+0.5))
        #initialize filter bank memory:
        initFB(N)
        #get and plot the synthesis impulse response of subband 0:
        y=np.ones((N,10))*0.0;
        y[0,0]=1.0
        xrek=np.ones(10*N)*0.0
        for m in range(10):
           xrek[(m*N):((m+1)*N)]=LDFBinv(y[:,m],fb)
        #Plots synthesis impulse response of subband 0:
        plt.plot(xrek)
        plt.title('MDCT Impulse Response of Subband 0 of the Synthesis FB')
        w,H=freqz(xrek,worN=2048)
        plt.figure()
        plt.plot(w,20*np.log10(np.abs(H)+1e-6))
        #Enlarge normalized frequencies in range of 0 to 0.1:
        plt.axis([0, 0.1, -60,5])
        plt.title('Its Enlarged Magnitude Frequency Response') 
        plt.ylabel('dB Attenuation')
        plt.xlabel('Normalized Frequency (pi is Nyquist Freq.)')
        plt.show()

        #initialize filter bank memory:
        initFB(N)

        CHUNK = N #Blocksize for the sound card
        WIDTH = 2 #2 bytes per sample
        CHANNELS = 1 #2
        RATE = 32000  #Sampling Rate in Hz
        p = pyaudio.PyAudio()

        a = p.get_device_count()
        print("device count=",a)

        for i in range(0, a):
            print("i = ",i)
            b = p.get_device_info_by_index(i)['maxInputChannels']
            print(b)
            b = p.get_device_info_by_index(i)['defaultSampleRate']
            print(b)
        #open sound device:
        stream = p.open(format=p.get_format_from_width(WIDTH),
                        channels=CHANNELS,
                        rate=RATE,
                        input=True,
                        output=True,
                        #input_device_index=3,
                        frames_per_buffer=CHUNK)


        print("* recording")
        #Waterfall diagram:
        #Size of waterfall diagramm:
        #N cols:
        rows=500
        cols=N
        
        frame=0.0*np.ones((rows,cols,3));
        ctr=0
        while(True):
            ctr=ctr+1
            #Reading from audio input stream into data with block length "CHUNK":
            data = stream.read(CHUNK)
            #Convert from stream of bytes to a list of short integers (2 bytes here) in "samples":
            #shorts = (struct.unpack( "128h", data ))
            shorts = (struct.unpack( 'h' * CHUNK, data ));
            samples=np.array(list(shorts),dtype=float);
            if np.max(np.abs(samples))>32000:
                print("Overload input")

            #shift "frame" 1 up:
            if (ctr%16 ==0):
               frame[0:(rows-1),:]=frame[1:rows,:]; 
           
            #This is the analysis MDCT of the input: 
            y=LDFB(samples[0:N],fb)

            #yfilt is the processed subbands, processing goes here:
            yfilt=y
            #yfilt=np.zeros(N)
            #yfilt[10:150]=y[10:150]
            #yfilt[30]=y[30]*4
            #yfilt[0:1024]=y[0:1024]
            
            #Waterfall color mapping when loop counter ctr reaches certain values:
            if ((display==True) & (ctr%16 ==0)):
               R=0.25*np.log((np.abs(yfilt/np.sqrt(N))+1))/np.log(10.0)
               #Red frame:
               frame[rows-1,:,2]=R
               #Green frame:
               frame[rows-1,:,1]=np.abs(1-2*R)
               #Blue frame:
               frame[rows-1,:,0]=1.0-R
               #frame[rows-1,:,0]=frame[rows-1,:,1]**3
               # Display the resulting frame
               cv2.imshow('MDCT Waterfall (end with "q")',frame)

            #Inverse, Synthesis filter bank:
            #Inverse/synthesis MDCT:
            xrek=LDFBinv(yfilt,fb);
            if np.max(np.abs(xrek))>32000:
                print("Overload output")
            xrek=np.clip(xrek, -32000, 32000)

            #Convert to short ints:
            xrek=np.array(xrek,dtype='int16')
            #converting from short integers to a stream of bytes in "data":
            #data=struct.pack('h' * len(samples), *samples);
            data=struct.pack('h' * len(xrek), *xrek);
            #Writing data back to audio output stream: 
            stream.write(data, CHUNK)

            #Keep window open until key 'q' is pressed:
            if cv2.waitKey(1) & 0xFF == ord('q'):
                break
        # When everything done, release the capture

        cv2.destroyAllWindows()

        stream.stop_stream()
        stream.close()
        p.terminate()

