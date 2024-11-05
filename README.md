
### guile-ffi-fftw

This is a minimal set of Guile FFI bindings for FFTW's ‘guru interface’,
  [www.fftw.org/fftw3_doc/Guru-Interface.html#Guru-Interface](https://www.fftw.org/fftw3_doc/Guru-Interface.html#Guru-Interface)

It provides two functions:

* `fftw-dft! transform-rank sign in out`
* `fftw-dft transform-rank sign in`

The array arguments `in` and `out` can be of type `c32` (complex float) or `c64` (complex double), and of any rank not smaller than `transform-rank`. The transform axes are last; please check the online help.

The bindings show how to interface between Guile and array libraries taking dense arrays as arguments. FFTW makes it easy by supporting arbitrary ranks and strides —which should be the norm!

The FFTW libraries are loaded with `dynamic-link`, either through the environment variable `GUILE_FFI_FFTW_LIBFFTW3_PATH`, or  on the default path. If only the float or only the double versions of the library are found, then `fftw-dft!` and `fftw-dft` will only support either `c32` or `c64` arguments.

These bindings being minimal, there is no support for computing & reusing plans, or split r/i transforms, or anything but straight complex DFTs. Contributions are welcome!

### Running the tests

The tests use SRFI-64.

```
$GUILE -L mod -s test/test-ffi-fftw.scm
```

### To do

- Arbitrary transform axes (make `transform-rank` into a list `transform-axes`).

### Links

- [https://www.fftw.org/](https://www.fftw.org/)
- [https://savannah.gnu.org/projects/guile-fftw](https://savannah.gnu.org/projects/guile-fftw) - set of bindings using a C extension
