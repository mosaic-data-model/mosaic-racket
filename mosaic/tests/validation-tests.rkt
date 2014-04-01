#lang racket

(require rackunit
         math/array
          "../validation.rkt"
         (prefix-in model: "../model.rkt"))

(provide validation-tests)

(define (m test)
  (cdr (first (messages test '()))))

(define (check-ok test)
  (check-pred empty? (messages test '())))

(define (check-message test expected-message)
  (let ([ms (messages test '())])
    (test-begin
     (check-equal? (length ms) 1)
     (check-equal? (cdr (first ms)) expected-message))))

; Fragment and universe  definitions used in the tests

(model:define-fragment water
                       #:atoms
                       (["H1" model:H]
                        ["H2" model:H]
                        ["O" model:O])
                       #:bonds
                       (["H1" "O" 'single]
                        ["H2" "O" 'single]))

(define water10
  (model:universe 'infinite '() "my-own"
                  (list (cons (model:fragment "water" water)  10))))

(define trunc-oct
  (model:universe 'cube
                  (list (cons (array #[1/2 1/2 1/2])
                              (array #[#[1 0 0] #[0 1 0] #[ 0 0 1]])))
                  "my-own"
                  (list (cons (model:fragment "water" water) 10))))

; Test suite

(define validation-tests
  (test-suite "validation"

; labels
(test-case "labels"

  (check-ok (validate-label "abc"))
  (check-ok (validate-label "ABC"))
  (check-ok (validate-label "[abc]+def-xyz*(012)/^"))

  (check-message (validate-label 42)
                 "label not a string: 42")

  (check-message (validate-label "")
                 "label is empty")

  (check-message (validate-label (make-string 40000 #\z))
                 "label too long: 40000 chars")

  (check-message (validate-label "a.b")
                 "label contains invalid character(s): \".\"")
  (check-message (validate-label "a b")
                 "label contains invalid character(s): \" \"")
  (check-message (validate-label "a;b:c")
                 "label contains invalid character(s): \";:\""))


; Atoms
(test-case "atoms"

  (check-ok (validate-atom (model:atom "C" model:C 1)))
  (check-ok (validate-atom (model:atom "H1"
                                       (model:atom-descr 'element "H")
                                       1)))
  (check-ok (validate-atom (model:atom "C"
                                       (model:atom-descr 'cgparticle "C-alpha")
                                       1)))
  (check-ok (validate-atom (model:atom "LP"
                                       (model:atom-descr 'dummy "lone-pair")
                                       1)))
  (check-ok (validate-atom (model:atom "a-foo"
                                       (model:atom-descr 'unknown "foo")
                                       1)))

  (check-message (validate-atom 42)
                 "not an atom")
  (check-message (validate-atom model:C)
                 "not an atom")

  (check-message (validate-atom (model:atom "C." model:C 1))
                 "label contains invalid character(s): \".\"")

  (check-message (validate-atom (model:atom "X"
                                            (model:atom-descr 'foo "foo")
                                            1))
                 "not a valid atom type: \"foo\"")

  (check-message (validate-atom (model:atom "H1"
                                            (model:atom-descr 'element "Xx")
                                            1))
                 "not a valid chemical element symbol: \"Xx\"")

  (check-message (validate-atom (model:atom "C" model:C 0))
                 "not a positive integer: 0")
  (check-message (validate-atom (model:atom "C" model:C 1.5))
                 "not a positive integer: 1.5"))

; Fragments
(test-case "fragments"

  (model:define-fragment err1.)

  (model:define-polymer foo 'polypeptide)

  (model:define-polymer poly-err 'polysomething)

  (model:define-fragment err2
                         #:atoms
                         (["H" model:H]
                          ["H" model:H]
                          ["O" model:O]))

  (model:define-fragment two-water
                         #:subfragments
                         (["mol1" water]
                          ["mol2" water])
                         #:bonds
                         (["mol1.O" "mol2.O" 'single]))

  (model:define-fragment two-water-err1
                         #:subfragments
                         (["mol1" water]
                          ["mol2" water])
                         #:bonds
                         (["mol1.O" "mol2" 'single]))

  (model:define-fragment two-water-err2
                         #:subfragments
                         (["mol1" water]
                          ["mol2" water])
                         #:bonds
                         (["mol1.O" "mol2.X" 'single]))

  (check-ok (validate-fragment (model:fragment "water" water)))
  (check-ok (validate-fragment (model:fragment "a-foo" foo)))
  (check-ok (validate-fragment (model:fragment "w2" two-water)))

  (check-message (validate-fragment 42)
                 "not a fragment")
  (check-message (validate-fragment water)
                 "not a fragment")

  (check-message (validate-fragment (model:fragment "water." water))
                 "label contains invalid character(s): \".\"")

  (check-message (validate-fragment (model:fragment "foo" err1.))
                 "species name contains invalid character(s): \".\"")

  (check-message (validate-fragment
                  (model:fragment "_" poly-err))
                 "not a valid polymer type: polysomething")

  (check-message (validate-fragment (model:fragment "water" err2))
                 "duplicate label(s)")

  (check-message (validate-fragment (model:fragment "w2" two-water-err1))
                 "node mol2 is not an atom")
  (check-message (validate-fragment (model:fragment "w2" two-water-err2))
                 "no atom mol2.X"))

; Universes
(test-case "universes"

  (check-ok (validate-universe water10))

  (check-ok (validate-universe trunc-oct))

  (check-message (validate-universe 42)
                 "not a universe")

  (check-message
   (validate-universe
    (model:universe 'some-shape '() "my-own"
                    (list (cons (model:fragment "water" water) 10))))
   "not a valid cell shape: \"some-shape\"")

  (check-message
   (validate-universe
    (model:universe 'infinite '() "a.b"
                    (list (cons (model:fragment "water" water) 10))))
   "convention contains invalid character(s): \".\"")

  (check-message
   (validate-universe
    (model:universe 'infinite '() "my-own"
                    (list (cons (model:fragment "water" water) 0))))
   "not a positive integer: 0")

  (check-message
   (validate-universe
    (model:universe 'infinite '() "my-own"
                    (list (cons 42 10))))
   "not a fragment")

  (check-message
   (validate-universe
    (model:universe 'infinite '() "my-own"
                    (list 'foo)))
   "not a fragment-count pair: foo"))

; Sequences
(test-case "sequences"

  (check-ok (validate-items (list water10 trunc-oct)))

  (check-message
   (validate-items (list 42))
   "unknown item type")

  (check-message
   (validate-items
    (list (model:universe 'infinite '() "my-own" (list 'foo))))
   "not a fragment-count pair: foo"))

; Configurations
(test-case "configurations"
  (check-ok (validate-configuration
             (model:configuration water10 (make-array #(30 3) 0.)
                                  (make-array #(0) (void)) 'float64)))

  (check-message
   (validate-configuration 42)
   "not a configuration")

  (check-message
   (validate-configuration
    (model:configuration 42 (make-array #(30 3) 0.)
                         (make-array #(0) (void)) 'float64))
   "not a universe: 42")

  (check-message
   (validate-configuration
    (model:configuration water10 (make-array #(30 3) 0.)
                         (make-array #(0) (void)) 'int8))
   "invalid number type: int8")

  (check-message
   (validate-configuration
    (model:configuration water10 42
                         (make-array #(0) (void)) 'float64))
   "not an array: 42")

  (check-message
   (validate-configuration
    (model:configuration water10 (make-array #(10 3) 0.)
                         (make-array #(0) (void)) 'float64))
   "invalid position shape #(10 3) (should be #(30 3))")

  (check-message
   (validate-configuration
    (model:configuration water10 (make-array #(30 3) 0.)
                         42 'float64))
   "not an array: 42")

  (check-message
   (validate-configuration
    (model:configuration water10 (make-array #(30 3) 0.)
                         (make-array #(1) 0.) 'float64))
   "invalid cell parameter shape #(1) (should be #(0))"))))
