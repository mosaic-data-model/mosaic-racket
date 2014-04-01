#lang racket/base
(require rackunit
         rackunit/text-ui
         "model-tests.rkt"
         "comparison-tests.rkt"
         "validation-tests.rkt"
         "xml-tests.rkt")

(provide all-tests
         main)

(define all-tests
  (make-test-suite "All MOSAIC tests"
    (list model-tests
          comparison-tests
          validation-tests
          xml-tests)))

(define (main) ;; run with "racket -tm <this-file>"
  (case 'text
    ((gui) ((dynamic-require 'rackunit/gui 'test/gui) all-tests #:wait? #t))
    ((text) (run-tests all-tests))))
