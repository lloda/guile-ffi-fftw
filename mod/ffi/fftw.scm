
; (c) Daniel Llorens - 2014, 2019

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU Lesser General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

;;; Commentary:
;; Access FFTW through Guile's FFI.
;;; Code:

(define-module (ffi fftw))
(import (system foreign) (srfi srfi-1))

; @TODO As an alternative go through installation.
(define libfftw3 (dynamic-link (let ((lpath (getenv "GUILE_FFI_FFTW_LIBFFTW3_PATH")))
                                 (if (and lpath (not (string=? lpath "")))
                                   (string-append lpath file-name-separator-string "libfftw3")
                                   "libfftw3"))))

;; http://www.fftw.org/doc/Guru-Complex-DFTs.html#Guru-Complex-DFTs
;; fftw_plan fftw_plan_guru_dft(int rank, const fftw_iodim *dims,
;;                              int howmany_rank, const fftw_iodim *howmany_dims,
;;                              fftw_complex *in, fftw_complex *out,
;;                              int sign, unsigned flags);

(define fftw_plan_guru_dft
  (pointer->procedure '* (dynamic-func "fftw_plan_guru_dft" libfftw3)
                      (list int '* int '* '* '* int unsigned-int)))
(define fftw_execute
  (pointer->procedure void (dynamic-func "fftw_execute" libfftw3)
                      (list '*)))
(define fftw_destroy_plan
  (pointer->procedure void (dynamic-func "fftw_destroy_plan" libfftw3)
                      (list '*)))

(define FFTW-ESTIMATE (ash 1 6))
(define FFTW-FORWARD -1)
(define FFTW-BACKWARD +1)

;; http://www.fftw.org/doc/Guru-vector-and-transform-sizes.html#Guru-vector-and-transform-sizes

(define make-iodims
  (let ((lim (ash 1 (+ -1 (* 8 (sizeof int))))))
    (lambda (in out)
      (append-map (lambda (s-in i-in s-out i-out)
                    (unless (= lim (max lim s-in i-in s-out i-out)) (throw 'fftw-sizes-too-large))
                    (unless (= s-in s-out) (throw 'fftw-mismatched-dimensions))
                    (list s-in i-in i-out))
                  (array-dimensions in)
                  (shared-array-increments in)
                  (array-dimensions out)
                  (shared-array-increments out)))))

(define (pick-iodims iodims n pick)
  (if (zero? n)
    (make-c-struct (list int) (list 0))
    (make-c-struct (make-list (* 3 n) int) (pick iodims (* 3 n)))))

(define (fftw-dft! k sign in out)
  "fftw-dft! k sign in out

Compute K-dimensional DFTs of array IN into array OUT, with given SIGN.  IN and
OUT must be 'c64 arrays of the same shape. For each K-cell of IN, a separate
K-DFT is computed and written in the corresponding cell of OUT.  SIGN is the
sign of the exponent of the transform and can be FFTW-FORWARD (-1) or
FFTW-BACKWARD (+1).  For example, if IN has shape (2 2 10 10), then after

 (fftw-dft! 2 FFTW-FORWARD IN OUT)

we have

 OUT[0, 0, ...] = 2D-DFT(IN[0, 0, ...])
 OUT[0, 1, ...] = 2D-DFT(IN[0, 1, ...])
 OUT[1, 0, ...] = 2D-DFT(IN[1, 0, ...])
 OUT[1, 1, ...] = 2D-DFT(IN[1, 1, ...])

unless OUT is the same array as IN, in which case OUT will contain the 2D-DFTs
of the original contents of IN. Otherwise it is assumed that OUT does not alias
IN.

See http://www.fftw.org/doc/FFTW-Reference.html for more information."

  (unless (<= 0 k (array-rank in))
    (throw 'fftw-bad-rank k (array-rank in)))
  (unless (= (array-rank out) (array-rank out))
    (throw 'fftw-mismatched-ranks (array-rank out) (array-rank out)))
  (unless (eq? 'c64 (array-type in) (array-type out))
    (throw 'fftw-bad-types (array-type in) (array-type out)))
  (unless (or (= +1 sign) (= -1 sign))
    (throw 'fftw-bad-sign sign))
  (let* ((iodims (make-iodims in out))
         (transform-k k)
         (repeat-k (- (array-rank in) k))
         (plan (fftw_plan_guru_dft
                transform-k (pick-iodims iodims transform-k take-right)
                repeat-k (pick-iodims iodims repeat-k take)
                (bytevector->pointer (shared-array-root in) (* (shared-array-offset in) (sizeof double) 2))
                (bytevector->pointer (shared-array-root out) (* (shared-array-offset out) (sizeof double) 2))
                sign FFTW-ESTIMATE)))
    (fftw_execute plan)
    (fftw_destroy_plan plan)))

(define (fftw-dft k sign in)
  "fftw-dft k sign in

Compute K-dimensional DFTs of 'c64 array IN, with given SIGN, and return the
result in a new 'c64 array of the same shape as IN.

See the documentation of fftw-dft! for more information."

  (let ((out (apply make-typed-array 'c64 *unspecified* (array-dimensions in))))
    (fftw-dft! k sign in out)
    out))

(export fftw-dft fftw-dft! FFTW-FORWARD FFTW-BACKWARD)
