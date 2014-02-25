
; Tests for (ffi fftw).
; (c) Daniel Llorens - 2014
; This is released and should depend only on (ffi fftw) and standard Guile.

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(import (ffi fftw) (srfi srfi-64) (srfi srfi-1))

; Test variable rank feature for (fftw-dft!).

(test-begin "fftw-dft! ones")
(define A (make-typed-array 'c64 1 2 2 2))
(define B1 (make-typed-array 'c64 1 2 2 2))
(fftw-dft! 1 +1 A B1)
(test-equal A (make-typed-array 'c64 1 2 2 2))
(test-equal B1 #3c64(((2 0) (2 0)) ((2 0) (2 0))))
(define B2 (make-typed-array 'c64 1 2 2 2))
(fftw-dft! 2 +1 A B2)
(test-equal A (make-typed-array 'c64 1 2 2 2))
(test-equal B2 #3c64(((4 0) (0 0)) ((4 0) (0 0))))
(define B3 (make-typed-array 'c64 1 2 2 2))
(fftw-dft! 3 +1 A B3)
(test-equal A (make-typed-array 'c64 1 2 2 2))
(test-equal B3 #3c64(((8 0) (0 0)) ((0 0) (0 0))))
(test-end "fftw-dft! ones")

; Test variable rank feature for (fftw-dft).

(test-begin "fftw-dft ones")
(define A (make-typed-array 'c64 1 2 2 2))
(test-equal (fftw-dft 1 +1 A) #3c64(((2 0) (2 0)) ((2 0) (2 0))))
(test-equal (fftw-dft 2 +1 A) #3c64(((4 0) (0 0)) ((4 0) (0 0))))
(test-equal (fftw-dft 3 +1 A) #3c64(((8 0) (0 0)) ((0 0) (0 0))))
(test-end "fftw-dft ones")

; Test variable rank feature for (fftw-dft) with nonstandard strides. FIXME:
; parameterize on input array.

(test-begin "fftw-dft ones, nonstandard strides I")
(define A (make-shared-array (make-typed-array 'c64 1.) (lambda i '()) 2 2 2))
(test-equal (fftw-dft 1 +1 A) #3c64(((2 0) (2 0)) ((2 0) (2 0))))
(test-equal (fftw-dft 2 +1 A) #3c64(((4 0) (0 0)) ((4 0) (0 0))))
(test-equal (fftw-dft 3 +1 A) #3c64(((8 0) (0 0)) ((0 0) (0 0))))
(test-end "fftw-dft ones, nonstandard strides I")

; Test signs.

(define i2pi (make-rectangular 0 (* 2 (acos -1))))

(define (delta n i)
  (let ((A (apply make-typed-array 'c64 0. n)))
    (apply array-set! A 1 i)
    A))

(define (delta-dft n i sign)
  (let ((A (apply make-typed-array 'c64 *unspecified* n)))
    (array-index-map!
     A (lambda k (exp (* sign i2pi (fold (lambda (i k n a) (+ a (* i k (/ n)))) 0 i k n)))))
    A))

(define (array-absolute-error a b)
  (and (equal? (array-dimensions a) (array-dimensions b))
       (let ((err 0))
         (array-for-each (lambda (a b) (set! err (max err (magnitude (- a b))))) a b)
         err)))

(define (delta-error sign . n)
  (let ((err 0.))
    (array-index-map!
     (apply make-shared-array (make-array #f) (lambda i '()) n)
     (lambda i
       (set! err (max err
                      (array-absolute-error (delta-dft n i sign)
                                            (fftw-dft (length n) sign (delta n i)))))))
    err))

(test-begin "fftw-dft deltas")
(test-approximate 0 (delta-error +1 4) 1e-15)
(test-approximate 0 (delta-error -1 4) 1e-15)
(test-approximate 0 (delta-error +1 3 4) 2e-15)
(test-approximate 0 (delta-error -1 3 4) 2e-15)
(test-approximate 0 (delta-error +1 2 3 4) 5e-15)
(test-approximate 0 (delta-error -1 2 3 4) 5e-15)
(test-end "fftw-dft deltas")
