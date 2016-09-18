This folder contains the files for [Peter Green's Nmix program](https://people.maths.bris.ac.uk/~mapjg/Nmix/) for univariate normal mixtures with a prior on the number of components. Nmix implements the reversible jump MCMC (RJMCMC) approach of Richardson and Green (1997), *Journal of the Royal Statistical Society: Series B*, 59(4), 731-792.

----------------------------------------------------------------------
## Contents

- Nmix.zip - Source files for Windows.
- Nmix.tar.gz - Source files for Unix-like environments (including Mac OS X and Linux).

----------------------------------------------------------------------
For Peter Green's instructions, see the readme.txt file in Nmix.zip or Nmix.tar.gz. 

Below is the process I used to compile Nmix on Mac and Windows.

### Nmix on Mac OS X

- If you haven't already done so, install the Command Line Tools (CLT) by opening a terminal window and entering: `xcode-select --install`.  In the window that pops up, click `Install` (**not** `Get Xcode`).
- Install [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS).
- Unpack [Nmix.tar.gz](Nmix.tar.gz) and run the following commands in the resulting folder:
```
gcc -c -o sd.o -DRETS -DSUNF sd.c
gfortran -O2 -o Nmix Nmix.f pnorm.f algama.f rgamma.f gauss4.f sd.o
```
You may get a message like `warning: implicit declaration of function 'time' is invalid in C99`. You can ignore this.
- Make sure Nmix is executable: `chmod +x Nmix`.

### Nmix on Windows

- If the pre-compiled executable works for you (`examples/Nmix.exe`), then you're good to go.  Otherwise, read on.
- Install gfortran by following the detailed instructions here: http://www.mingw.org/wiki/Getting_Started. Choose to install mingw32-base, mingw32-gcc-fortran, and msys-base.
- Unzip Nmix.zip and copy the contents to your MSYS home directory, e.g., mine is `C:\MinGW\msys\1.0\home\jeff\`.
- Start MSYS and run the following commands in the folder containing the Nmix source files:
```
gcc -c -o sd.o -DRETS -DSUNF sd.c
gfortran -O2 -static -o Nmix Nmix.f pnorm.f algama.f rgamma.f gauss4.f sd.o
```
You may get a message like `warning: implicit declaration of function 'time' is invalid in C99`. You can ignore this.



