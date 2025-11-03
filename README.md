# FFTOD
Fftod: is a packge for computing the Fast Fourier Transform (FFT) written in the Odin language.
It has been created primarily to serve as the internal fft engine for [Nuod](https://github.com/Omair-R/Nuod). It doesn't try to be the fastest, but it does provide decent performance with a convenient interface. 

## Features
- Operates with signals of any length by utilizing mixed-radix [Stockham](https://ieeexplore.ieee.org/document/1447887) (2, 3, 5) and [Bluestein](https://ieeexplore.ieee.org/document/1162132) implementations.
- It has no dependances outside of the core and intrinsic Odin libraries[^1].
- It provides several interfaces with varying levels of control and convenience.
- Works with all the complex data types provided by Odin (through parametric polymorphism).
- In-place (in-situ) and out-of-place (ex-situ) variants.

[^1]: the "tests" folder has a dependency on the Ooura library, purely to benchmark against. Nevertheless, Fftod itself is free of external dependencies and can be dropped into your project without worry. 

## Usage
Adding the library to your project should work by copying the "fftod" to your project repository. 

To use the library for performing the FFT, the following examples have been prepared

> [!NOTE]
> the variable "x" is a slice of complex values with a size of N that represents your signal. For these examples it is replaces with a random slice as per:
```odin
x := make([]complex128, N)		
defer delete(x)

for i in 0..<N{
  x[i] = complex(rand.float64(), rand.float64())
}
```
### High-level Interface
```odin
// perform the fft out-of-place (allocates the array f)
f := fftod.out_fft(x, allocator=..., location=...) // allocator and location are optional.
defer delete(f)
// perform the inverse fft out-of-place
x_ := fftod.out_ifft(f, allocator=..., location=...)
defer delete(x_)

// perform the fft in-place (The allocator is required for work buffers)
fftod.in_fft(x, allocator=..., location=...)

// perform the inverse fft in-place
fftod.in_ifft(x, allocator=..., location=...)
```

### Lower-level Interface
```odin
// perform the fft out-of-place with the specific algorithm suitable for the length of signal (N)
f := fftod.out_fft_stockham_radix4(x, inverse=false, allocator=..., location=...) // assuming the input is a power of 4.
defer delete(f)
// perform the inverse fft out-of-place
x_ := fftod.out_fft_stockham_radix4(f, inverse=true, allocator=..., location=...)  
defer delete(x_)

// perform the fft in-place with the specific algorithm suitable for the length of signal (N)
fftod.in_fft_stockham_radix3(x, inverse=false, allocator=..., location=...) // assuming the input is a power of 3.
// perform the inverse fft in-place
fftod.in_fft_stockham_radix3(x, inverse=true, allocator=..., location=...)
```

### Lowest-level Interface
```odin
// Create a plan instance that maintains all work buffers for reuse with several signals of the same length. 
plan := fftod.make_bluestein_plan(N, complex128, allocator:=..., location:=...)
defer delete_bluestein_plan(plan)

// you may reuse "plan" as many times as you which, so long as the length of each signal is equal to N.
f := fftod.out_fft_bluestein(x, plan, inverse=false, allocator=..., location=...)
defer delete(f)
```
## Additional Note
It might seem a bit odd to use the Stockham algorithm instead of the common Cooley–Tukey; However, the former seemed to deliver a better performance all while providing a nicer iterative form for building the mixed-radix version. It also does not require any digit or bit reversal.

## References 
- [Postpischil, E. (2004). Construction of a High-Performance FFT. *Mathematics, Design, and Implementation Guide*](https://edp.org/work/Construction.pdf)
- [Brigham, E. (1988). Fast Fourier Transform and its applications.](https://espy.folk.ntnu.no/Second/Master_Student_Resources/reference_books/Bingham_The_FFT/FFT%20and%20its%20applications.pdf)
- [Lyons, R. G. (1996). Understanding digital signal processing.](https://www.mikrocontroller.net/attachment/341426/Understanding_digital_signal_processing.pdf)
- [Bluestein, L.  (1970). A linear filtering approach to the computation of discrete Fourier transform *IEEE Transactions on Audio and Electroacoustics*](https://ieeexplore.ieee.org/document/1162132)
- [Cochran, W. T. (1967). What is the fast Fourier transform? *Proceedings of the IEEE*](https://ieeexplore.ieee.org/document/1447887)
- [Introduction to the Stockham FFT *OTFFT: High Speed FFT library*](http://wwwa.pikara.ne.jp/okojisan/otfft-en/stockham1.html)
- [Burrus, C. S. The Chirp Z-Transform or Bluestein's Algorithm](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Fast_Fourier_Transforms_(Burrus)/04%3A_The_DFT_as_Convolution_or_Filtering/4.03%3A_The_Chirp_Z-Transform_or_Bluestein's_Algorithm)
### Original License for The OOURA Package
```
Copyright(C) 1996-2001 Takuya OOURA
email: ooura@mmm.t.u-tokyo.ac.jp
download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
You may use, copy, modify this code for any purpose and 
without fee. You may distribute this ORIGINAL package.
```
