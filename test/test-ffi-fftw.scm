; -*- mode: scheme; coding: utf-8 -*-

; (c) Daniel Llorens - 2014-2024
; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU Lesser General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

;;; Commentary:
;; Tests for (ffi fftw).
;;; Code:

(import (ffi fftw) (srfi srfi-64) (srfi srfi-1) (ice-9 match))

(test-begin "ffi-fftw")

; Various sorts of arrays.

(define (make-A-compact type)
  (make-typed-array type 0 2 2 2))

(define (make-A-strided type)
  (make-shared-array (make-typed-array type 0 8 8 8) (lambda i (map * i '(2 3 4))) 2 2 2))

(define (make-A-offset type)
  (make-shared-array (make-typed-array type 0 8 8 8) (lambda i (map + i '(2 3 4))) 2 2 2))

; Test variable rank feature for (fftw-dft!).

(define (case-fftw-dft!-ones type tag make-A make-B)
  (let ((case-name (format #f "fftw-dft! ones, ~a" tag))
        (ref (make-typed-array type 1 2 2 2))
        (A (make-A type))
        (B (make-B type)))
    (array-fill! A 1.)
    (test-begin case-name)
    (array-fill! B 11.)
    (fftw-dft! 1 +1 A B)
    (test-equal B (list->typed-array type 3 '(((2 0) (2 0)) ((2 0) (2 0)))))
    (test-equal A ref)
    (array-fill! B 22.)
    (fftw-dft! 2 +1 A B)
    (test-equal B (list->typed-array type 3 '(((4 0) (0 0)) ((4 0) (0 0)))))
    (test-equal A ref)
    (array-fill! B 33.)
    (let ((C (fftw-dft! 3 +1 A B)))
      (test-equal B (list->typed-array type 3 '(((8 0) (0 0)) ((0 0) (0 0)))))
      (test-equal A ref)
      (test-eq B C))
    (test-end case-name)))

; Test variable rank feature for (fftw-dft).

(define (case-fftw-dft-ones type tag make-A)
  (let ((case-name (format #f "fftw-dft ones, ~a" tag))
        (A (make-A type)))
    (array-fill! A 1.)
    (test-begin case-name)
    (test-equal (fftw-dft 1 +1 A) (list->typed-array type 3 '(((2 0) (2 0)) ((2 0) (2 0)))))
    (test-equal (fftw-dft 2 +1 A) (list->typed-array type 3 '(((4 0) (0 0)) ((4 0) (0 0)))))
    (test-equal (fftw-dft 3 +1 A) (list->typed-array type 3 '(((8 0) (0 0)) ((0 0) (0 0)))))
    (test-end case-name)))

(for-each
    (lambda (type)
      (when (match type ('c32 libfftw3f) ('c64 libfftw3))
        (test-begin (symbol->string type))
        (case-fftw-dft!-ones type "compact-compact" make-A-compact make-A-compact)
        (case-fftw-dft!-ones type "compact-strided" make-A-compact make-A-strided)
        (case-fftw-dft!-ones type "compact-offset" make-A-compact make-A-offset)
        (case-fftw-dft!-ones type "strided-compact" make-A-strided make-A-compact)
        (case-fftw-dft!-ones type "strided-strided" make-A-strided make-A-strided)
        (case-fftw-dft!-ones type "strided-offset" make-A-strided make-A-offset)
        (case-fftw-dft!-ones type "offset-compact" make-A-offset make-A-compact)
        (case-fftw-dft!-ones type "offset-strided" make-A-offset make-A-strided)
        (case-fftw-dft!-ones type "offset-offset" make-A-offset make-A-offset)

        (case-fftw-dft-ones type "fresh array" make-A-compact)
        (case-fftw-dft-ones type "strided" make-A-strided)
        (case-fftw-dft-ones type "offset" make-A-offset)
        (case-fftw-dft-ones type "reshaped scalar" (lambda (type) (make-shared-array (make-typed-array type 0.) (lambda i '()) 2 2 2)))
        (test-end (symbol->string type))))
  '(c32 c64))

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
