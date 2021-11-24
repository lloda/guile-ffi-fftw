
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
  "Compute @var{k}-dimensional DFTs of array @var{in} with given @var{sign}, and
write them to array @var{out}.  @var{in} and @var{out} must be c64 arrays of the
same shape. For each @var{k}-cell of @var{in}, a separate @var{k}-DFT is
computed and written to the matching cell of @var{out}.  @var{sign} is the sign
of the exponent of the transform and can be @code{FFTW-FORWARD} (-1) or
@code{FFTW-BACKWARD} (+1).  For example, if @var{in} has shape @code{(2 2 10
10)}, then after

@lisp
 (fftw-dft! 2 FFTW-FORWARD IN OUT)
@end lisp

we have

@example
 out[0, 0, ...] = 2D-DFT(in[0, 0, ...])
 out[0, 1, ...] = 2D-DFT(in[0, 1, ...])
 out[1, 0, ...] = 2D-DFT(in[1, 0, ...])
 out[1, 1, ...] = 2D-DFT(in[1, 1, ...])
@end example

If @var{out} is the same array as @var{in}, then @var{out} will contain the
2D-DFTs of the original contents of @var{in}. Otherwise it is assumed that
@var{out} does not alias @var{in} and if it does, then the contents of both
@var{out} and @var{in} become unspecified.

This function returns the output array @var{out}.

See http://www.fftw.org/doc/FFTW-Reference.html for more information.

See also: fttw-dft
"

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
    (fftw_destroy_plan plan)
    out))

(define (fftw-dft k sign in)
  "Compute @var{k}-dimensional DFTs of 'c64 array @var{in}, with given
@var{sign}, and return the result in a new 'c64 array of the same shape as
@var{in}.

See also: fttw-dft!
"

  (let ((out (apply make-typed-array 'c64 *unspecified* (array-dimensions in))))
    (fftw-dft! k sign in out)))

(export fftw-dft fftw-dft! FFTW-FORWARD FFTW-BACKWARD)
