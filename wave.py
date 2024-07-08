import math
import matplotlib.pyplot as plt
import numpy as np
from mpmath import *

fs = 1000000.0
f0 = 1000.0
phi0 = 0.0

PERIODS = 27

omega = 2.0 * f0 * math.pi
k = 155.00
kmin = 10.0
kstep = 1.0

kminarg = 1.0
kminval = 1e9

times = np.empty(round(PERIODS * fs / f0))
wavre = np.empty(round(PERIODS * fs / f0), np.float64)
wavim = np.empty(round(PERIODS * fs / f0), np.float64)

wavk = np.empty(100)
wavsum = np.empty(100)

mp.prec = 128

def fillwavelet(k):
    for i in range(0, (round(PERIODS * fs / f0))):
        t: float
        ti = i - (round(PERIODS / 2 * fs / f0))
        times[i] = ti
        t = float(ti) * omega / fs
        wavre[i] = math.cos(t) * math.exp(-t * t / k)
        wavim[i] = math.sin(t) * math.exp(-t * t / k)


def findk(kmin):
    k = kmin
    global kminval
    global kminarg
    fillwavelet(k)
    kminval = sum(wavre)
    for i in range(0, 100):
        fillwavelet(k)
        wavk[i] = k
        wavsum[i] = sum(wavre)
        if abs(wavsum[i]) < abs(kminval):
            kminval = abs(wavsum[i])
            kminarg = k
        k = k + kstep


fillwavelet(kmin)

Sum = sum(wavre)
print('wavre:', Sum)

# plt.plot(times, wavre)
# plt.plot(times, wavim)
# plt.show()

findk(kmin)
for i in range(5):
    kmin = kminarg - kstep
    kstep = kstep * 0.1
    findk(kmin)
    print('Iteration ', i + 1, ' : ', sum(wavre))

fillwavelet(kminarg)


def test_phase():
    delta_sqr = 0.0
    for i in range(90):
        phi = i * math.pi / 180
        wav_re: float = 0
        wav_im: float = 0
        pw_re: int = 0
        pw_im: int = 0
        for j in range(0, (round(PERIODS * fs / f0))):
            x = round((2**16 - 1) * math.cos(2.0 * math.pi * (j - (round(PERIODS / 2 * fs / f0))) * f0 / fs - phi))
            wav_re += wavre[j] * x
            wav_im += wavim[j] * x
            pw_re += wavre[j]**2
            pw_im += wavim[j]**2
        phase = math.atan(wav_im / wav_re)
        delta = phase - phi
        delta_sqr += sqrt(delta * delta / 90)
        print(i, delta)
    print("Sum of delta sqr = ", delta_sqr)

def test_phase_fourier():
    wampl = 2**32
    PERIODS = 27
    delta_sqr = 0.0
    for i in range(90):
        phi: float = math.radians(i)
        wav_re: int = 0
        wav_im: int = 0
        pw_re: int = 0
        pw_im: int = 0
        for j in range(0, (round(PERIODS * fs / f0))):
            x = round((2**16 - 1) * mp.cos(2.0 * math.pi * (j - (round(PERIODS / 2 * fs / f0))) * f0 / fs - phi))
            wav_re += x * round(wampl * mp.cos(2.0 * math.pi * (j - (round(PERIODS / 2 * fs / f0))) * f0 / fs))
            wav_im += x * round(wampl * mp.sin(2.0 * math.pi * (j - (round(PERIODS / 2 * fs / f0))) * f0 / fs))
            pw_re += round(wampl * np.cos(2.0 * math.pi * (j - (round(PERIODS / 2 * fs / f0))) * f0 / fs))**2
            pw_im += round(wampl * np.sin(2.0 * math.pi * (j - (round(PERIODS / 2 * fs / f0))) * f0 / fs))**2
        phase = math.atan(wav_im / wav_re)
        delta = phase - phi
        delta_sqr += sqrt(delta * delta / 90)
        print(i, wav_re, wav_im, delta)
    print('Power of kernel: ', pw_re, pw_im)
    print("Sum of delta sqr = ", delta_sqr)


# print(kminarg, kminval)
# plt.plot(wavk, wavsum)
# plt.show()

test_phase()
test_phase_fourier()
