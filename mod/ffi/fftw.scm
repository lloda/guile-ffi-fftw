; -*- mode: scheme; coding: utf-8 -*-

; (c) Daniel Llorens - 2014-2024
; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU Lesser General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

;;; Commentary:
;; Access FFTW through Guile's FFI.
;;; Code:

(define-module (ffi fftw)
  #:export (libfftw3 libfftw3f fftw-dft fftw-dft! FFTW-FORWARD FFTW-BACKWARD))
(import (system foreign) (srfi srfi-1) (srfi srfi-71) (ice-9 match))

; @TODO As an alternative go through installation.
(define (dylink name)
  (catch 'misc-error
    (lambda ()
      (dynamic-link (let ((lpath (getenv "GUILE_FFI_FFTW_LIBFFTW3_PATH")))
                      (if (and lpath (not (string=? lpath "")))
                        (string-append lpath file-name-separator-string name)
                        name))))
    (lambda (k . args)
      (format (current-error-port) "Can't find ~a [~a ~a]\n" name k args)
      #f)))

(define libfftw3 (dylink "libfftw3"))
(define libfftw3f (dylink "libfftw3f"))
(unless (or libfftw3 libfftw3f)
  (throw 'fftw-libraries-cannot-be-found))

;; http://www.fftw.org/doc/Guru-Complex-DFTs.html#Guru-Complex-DFTs
;; fftw_plan fftw_plan_guru_dft(int rank, const fftw_iodim *dims,
;;                              int howmany_rank, const fftw_iodim *howmany_dims,
;;                              fftw_complex *in, fftw_complex *out,
;;                              int sign, unsigned flags);

(define-values (fftw_plan_guru_dft fftw_execute fftw_destroy_plan)
  (if libfftw3
    (values
       (pointer->procedure '* (dynamic-func "fftw_plan_guru_dft" libfftw3)
                           (list int '* int '* '* '* int unsigned-int))
       (pointer->procedure void (dynamic-func "fftw_execute" libfftw3)
                      (list '*))
       (pointer->procedure void (dynamic-func "fftw_destroy_plan" libfftw3)
                           (list '*)))
    (values #f #f #f)))

(define-values (fftwf_plan_guru_dft fftwf_execute fftwf_destroy_plan)
  (if libfftw3f
    (values
     (pointer->procedure '* (dynamic-func "fftwf_plan_guru_dft" libfftw3f)
                         (list int '* int '* '* '* int unsigned-int))
     (pointer->procedure void (dynamic-func "fftwf_execute" libfftw3f)
                         (list '*))
     (pointer->procedure void (dynamic-func "fftwf_destroy_plan" libfftw3f)
                         (list '*)))
    (values #f #f #f)))

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
write them to array @var{out}. @var{in} and @var{out} must be arrays of the
same shape and type, either 'c32 or 'c64. For each @var{k}-cell of @var{in}, a
separate @var{k}-DFT is computed and written to the matching cell of @var{out}.
@var{sign} is the sign of the exponent of the transform and can be
@code{FFTW-FORWARD} (-1) or @code{FFTW-BACKWARD} (+1). For example, if @var{in}
has shape @code{(2 2 10 10)}, then after

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
2D-DFTs of the original contents of @var{in}. If @var{out} aliases @var{in}
without being the same, the contents of both @var{out} and @var{in} become
undefined.

This function returns the output array @var{out}.

See http://www.fftw.org/doc/FFTW-Reference.html for more information.

See also: fttw-dft
"
  (define type (array-type in))
  (define rank (array-rank in))
  (unless (and (eq? type (array-type out)) (or (eq? type 'c64) (eq? type 'c32)))
    (throw 'fftw-bad-types type (array-type out)))
  (unless (<= 0 k rank)
    (throw 'fftw-bad-rank k rank))
  (unless (= rank (array-rank out))
    (throw 'fftw-mismatched-ranks rank (array-rank out)))
  (unless (or (= +1 sign) (= -1 sign))
    (throw 'fftw-bad-sign sign))

  (let* ((iodims (make-iodims in out))
         (transform-k k)
         (repeat-k (- rank k))
         (real plan execute destroy
           (match type
             ('c32 (unless libfftw3f (throw 'fftw-is-not-available-for-c32))
                   (values float fftwf_plan_guru_dft fftwf_execute fftwf_destroy_plan))
             ('c64 (unless libfftw3 (throw 'fftw-is-not-available-for-c64))
                   (values double fftw_plan_guru_dft fftw_execute fftw_destroy_plan))))
         (p (plan
             transform-k (pick-iodims iodims transform-k take-right)
             repeat-k (pick-iodims iodims repeat-k take)
             (bytevector->pointer (shared-array-root in) (* (shared-array-offset in) (sizeof real) 2))
             (bytevector->pointer (shared-array-root out) (* (shared-array-offset out) (sizeof real) 2))
             sign FFTW-ESTIMATE)))
    (execute p)
    (destroy p)
    out))

(define (fftw-dft k sign in)
  "Compute @var{k}-dimensional DFTs of 'c64 or 'c32 array @var{in}, with given
@var{sign}, and return the result in a new array of the same shape and type as
@var{in}.

See also: fttw-dft!
"
  (let ((out (apply make-typed-array (array-type in) *unspecified* (array-dimensions in))))
    (fftw-dft! k sign in out)))
