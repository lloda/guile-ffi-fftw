
; Tests for (ffi fftw).
; (c) Daniel Llorens - 2014
; This is released and should depend only on (ffi fftw) and standard Guile.

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU Lesser General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(import (ffi fftw) (srfi srfi-64) (srfi srfi-1))

(set! test-log-to-file #f)
(test-begin "ffi-fftw")

; Various sorts of arrays.

(define (make-A-compact)
  (make-typed-array 'c64 0 2 2 2))

(define (make-A-strided)
  (make-shared-array (make-typed-array 'c64 0 8 8 8) (lambda i (map * i '(2 3 4))) 2 2 2))

(define (make-A-offset)
  (make-shared-array (make-typed-array 'c64 0 8 8 8) (lambda i (map + i '(2 3 4))) 2 2 2))

; Test variable rank feature for (fftw-dft!).

(define (case-fftw-dft!-ones tag make-A make-B)
  (let ((case-name (format #f "fftw-dft! ones, ~a" tag))
        (ref (make-typed-array 'c64 1 2 2 2))
        (A (make-A))
        (B (make-B)))
    (array-fill! A 1.)
    (test-begin case-name)
    (array-fill! B 11.)
    (fftw-dft! 1 +1 A B)
    (test-equal B #3c64(((2 0) (2 0)) ((2 0) (2 0))))
    (test-equal A ref)
    (array-fill! B 22.)
    (fftw-dft! 2 +1 A B)
    (test-equal B #3c64(((4 0) (0 0)) ((4 0) (0 0))))
    (test-equal A ref)
    (array-fill! B 33.)
    (fftw-dft! 3 +1 A B)
    (test-equal B #3c64(((8 0) (0 0)) ((0 0) (0 0))))
    (test-equal A ref)
    (test-end case-name)))

(case-fftw-dft!-ones "compact-compact" make-A-compact make-A-compact)
(case-fftw-dft!-ones "compact-strided" make-A-compact make-A-strided)
(case-fftw-dft!-ones "compact-offset" make-A-compact make-A-offset)
(case-fftw-dft!-ones "strided-compact" make-A-strided make-A-compact)
(case-fftw-dft!-ones "strided-strided" make-A-strided make-A-strided)
(case-fftw-dft!-ones "strided-offset" make-A-strided make-A-offset)
(case-fftw-dft!-ones "offset-compact" make-A-offset make-A-compact)
(case-fftw-dft!-ones "offset-strided" make-A-offset make-A-strided)
(case-fftw-dft!-ones "offset-offset" make-A-offset make-A-offset)

; Test variable rank feature for (fftw-dft).

(define (case-fftw-dft-ones tag A)
  (let ((case-name (format #f "fftw-dft ones, ~a" tag)))
    (array-fill! A 1.)
    (test-begin case-name)
    (test-equal (fftw-dft 1 +1 A) #3c64(((2 0) (2 0)) ((2 0) (2 0))))
    (test-equal (fftw-dft 2 +1 A) #3c64(((4 0) (0 0)) ((4 0) (0 0))))
    (test-equal (fftw-dft 3 +1 A) #3c64(((8 0) (0 0)) ((0 0) (0 0))))
    (test-end case-name)))

(case-fftw-dft-ones "fresh array" (make-A-compact))
(case-fftw-dft-ones "strided" (make-A-strided))
(case-fftw-dft-ones "offset" (make-A-offset))
(case-fftw-dft-ones "reshaped scalar" (make-shared-array (make-typed-array 'c64 0.) (lambda i '()) 2 2 2))

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

(define error-count (test-runner-fail-count (test-runner-current)))
(test-end "ffi-fftw")
(exit error-count)
