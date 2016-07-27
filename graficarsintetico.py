
import matplotlib.pyplot as plt
import scipy.signal as signal
# import matplotlib.patches as ptch
import numpy as np
from matplotlib import rcParams
import matplotlib
import sys
import sounddevice as sd
import time
fn=44150        
a1, a2 = signal.butter(2, 5/fn, btype='highpass', analog=False, output='ba')
###Sonograma y vs A B DIA
vs=np.loadtxt(sys.argv[1])
sintetizado=np.loadtxt('sintetizado.'+sys.argv[1]+'.dat')
envolvente=np.loadtxt('envolvente.'+sys.argv[1]+'.dat')
filtrado1=np.loadtxt('ptraqout.'+sys.argv[1]+'.dat')
pfinal=np.loadtxt('pfinal.'+sys.argv[1]+'.dat')
envcanto=np.loadtxt('envolvente.'+sys.argv[2]+'.dat')
canto=np.loadtxt(sys.argv[2])
noche=np.loadtxt('pfinal.'+sys.argv[3]+'.dat')
# sintetizado=signal.filtfilt(a1,a2,sintetizado)
filtrado1=signal.filtfilt(a1,a2,filtrado1)
a1, a2 = signal.butter(2, 400/fn, btype='highpass', analog=False, output='ba')
pfinal=signal.filtfilt(a1,a2,pfinal)
noche=signal.filtfilt(a1,a2,noche)

noche/=np.max(np.abs(noche))
sintetizado/=np.max(np.abs(sintetizado))
filtrado1/=np.max(np.abs(filtrado1))
pfinal/=np.max(np.abs(pfinal))


Fs=44150
t=np.linspace(0,len(sintetizado)/Fs,len(sintetizado))

print("sintetizado max=",np.max(sintetizado))
print("sintetizado min=",np.min(sintetizado))

if(np.max(sintetizado)==np.min(sintetizado)):
    sys.exit("MINIMO=MAXIMO")

plt.figure()
plt.title('envolvente normalizada')
plt.plot(t,envolvente)

plt.figure()
plt.title('envolvente canto')
plt.plot(t,envcanto)


plt.figure(figsize=(20,10),dpi=80)
plt.plot(t[0:len(vs)],vs)


f, (ax1, ax2) = plt.subplots(2, sharex=True,figsize=(20,10),dpi=80)
ax1.plot(t[0:len(vs)],vs)
ax1.set_title('vs')
# ax2.specgram(sintetizado,pad_to=600,vmin=vmin,vmax=vmax,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)
ax2.plot(t,sintetizado)
ax2.set_ylim([np.min(sintetizado),np.max(sintetizado)])


cmap = plt.get_cmap('Greys')
# vmin = 20*np.log10(np.max(sintetizado))-100
# vmax=-40

# rcParams.update({'font.size': 32})
# fig=plt.figure(num=None, figsize=(20,10), dpi=80,  edgecolor='k')
# axes = plt.gca()
# Pxx, freqs, bins, im = plt.specgram(sintetizado,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)
# axes.set_ylim([0,10000])
# plt.xlabel("Tiempo (s)")

# scale = 1e3                     # KHz
# ticks = matplotlib.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale))
# axes.yaxis.set_major_formatter(ticks)
# plt.ylabel("Frecuencia (kHz)")


f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(20,10))
ax1.plot(t[0:len(vs)],vs)
ax1.set_title('vs')
ax2.specgram(sintetizado,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)
ax2.set_ylim([0,10000])


# sd.play(sintetizado, 44150)

fig=plt.figure(num=None, figsize=(20,10), dpi=80,  edgecolor='k')
plt.title("FILTRADO1")
axes = plt.gca()
Pxx, freqs, bins, im = plt.specgram(filtrado1,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)
axes.set_ylim([0,10000])
plt.xlabel("Tiempo (s)")
time.sleep(1)
sd.play(filtrado1, 44150)


f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(20,10))
plt.suptitle('Pfinal')
ax1.plot(t[0:len(pfinal)],pfinal)

Pxx, freqs, bins, im = ax2.specgram(pfinal,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)
axes = plt.gca()
axes.set_ylim([0,10000])
plt.xlabel("Tiempo (s)")



sd.play(canto,44150)
time.sleep(1)

sd.play(pfinal, 44150)
time.sleep(1)

sd.play(noche, 44150)

f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(20,10))
plt.suptitle('SONOGRAMAS')
Pxx, freqs, bins, im = ax1.specgram(canto,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)
ax1.set_ylim([0,10000])
Pxx, freqs, bins, im = ax2.specgram(pfinal,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)

ax2.set_ylim([0,10000])
plt.xlabel("Tiempo (s)")

Pxx, freqs, bins, im = ax3.specgram(noche,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)

ax3.set_ylim([0,10000])
plt.show(block=True)
