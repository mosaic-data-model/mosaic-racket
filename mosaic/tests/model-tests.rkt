#lang racket

(require rackunit
         math/array
         "../model.rkt"
         (prefix-in interface: "../interface.rkt"))

(provide model-tests)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Universes
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-fragment water
  #:atoms
  (["H1" H]
   ["H2" H]
   ["O" O])
  #:bonds
  (["H1" "O" 'single]
   ["H2" "O" 'single]))

(define-fragment alanine-nter
  #:atoms
  (["H1" H]
   ["H2" H]
   ["H3" H]
   ["N" N]
   ["CA" C]
   ["HA" H]
   ["C" C]
   ["O" O]
   ["CB" C]
   ["HB1" H]
   ["HB2" H]
   ["HB3" H])
  #:bonds
  (["N" "H1" 'single]
   ["N" "H2" 'single]
   ["N" "H3" 'single]
   ["N" "CA" 'single]
   ["CA" "HA" 'single]
   ["CA" "C" 'single]
   ["C" "O" "double"]
   ["CA" "CB" 'single]
   ["CB" "HB1" 'single]
   ["CB" "HB2" 'single]
   ["CB" "HB3" 'single]))

(define-fragment alanine-cter
  #:atoms
  (["H" H]
   ["N" N]
   ["CA" C]
   ["HA" H]
   ["C" C]
   ["O" O]
   ["OXT" O]
   ["CB" C]
   ["HB1" H]
   ["HB2" H]
   ["HB3" H])
  #:bonds
  (["N" "H" 'single]
   ["N" "CA" 'single]
   ["CA" "HA" 'single]
   ["CA" "C" 'single]
   ["C" "O" 'single]
   ["C" "OXT" 'single]
   ["CA" "CB" 'single]
   ["CB" "HB1" 'single]
   ["CB" "HB2" 'single]
   ["CB" "HB3" 'single]))

(define-polymer di-alanine 'polypeptide
  #:subfragments
  (["ALA1" alanine-nter]
   ["ALA2" alanine-cter])
  #:bonds
  (["ALA1.C" "ALA2.N" 'single]))

(define water10 (universe 'infinite '() "my-own"
                          (list (cons (fragment "water" water) 10))))

(define trunc-oct (universe 'cube
                            (list (cons (array #[1/2 1/2 1/2])
                                        (array #[#[1 0 0] #[0 1 0] #[ 0 0 1]])))
                            "my-own"
                            (list (cons (fragment "water" water) 10))))


(define model-tests
  (test-suite "model"

; Test Mosaic interface

(test-case "interface"

  (let ([a (dict-ref water "O")])
    (check-equal? (interface:node.label a) "O")
    (check-equal? (interface:atom.type a) "element")
    (check-equal? (interface:atom.name a) "O") 
    (check-equal? (interface:atom.nsites a) 1))

  (let ([w (fragment "w" water)])
    (check-equal? (interface:node.label w) "w")
    (check-equal? (interface:fragment.species w) "water")
    (check-equal? (interface:fragment.subfragments w) '())
    (check-equal? (interface:fragment.atoms w)
                  (list (atom "H1" H 1) (atom "H2" H 1) (atom "O" O 1)))
    (check-equal? (interface:fragment.lookup-node w "O") (atom "O" O 1))
    (check-equal? (interface:fragment.lookup-node w "H1") (atom "H1" H 1))
    (check-equal? (interface:fragment.lookup-node w "H2") (atom "H2" H 1))
    (check-equal? (interface:fragment.bonds w)
                  (list #("H1" "O" "single") #("H2" "O" "single")))
    (check-equal? (interface:fragment.polymer? w) #f)
    (check-exn exn:fail:contract? (λ () (interface:fragment.polymer-type w))))

  (let ([da (fragment "da" di-alanine)])
    (check-equal? (interface:node.label da) "da")
    (check-equal? (interface:fragment.species da) "di-alanine")
    (check-equal? (interface:fragment.subfragments da)
                  (list (fragment "ALA1" alanine-nter)
                        (fragment "ALA2" alanine-cter)))
    (check-equal? (interface:fragment.atoms da) '())
    (check-equal? (interface:fragment.lookup-node da "ALA1")
                  (fragment "ALA1" alanine-nter))
    (check-equal? (interface:fragment.bonds da)
                  (list #("ALA1.C" "ALA2.N" "single")))
    (check-equal? (interface:fragment.polymer? da) #t)
    (check-equal? (interface:fragment.polymer-type da) "polypeptide"))

  (check-equal? (interface:universe.cell-shape water10) "infinite")
  (check-equal? (interface:universe.symmetry-transformations water10) '())
  (check-equal? (interface:universe.convention water10) "my-own")
  (check-equal? (interface:universe.molecules water10)
                (list (cons (fragment "water" water) 10)))

  (check-equal? (interface:universe.cell-shape trunc-oct) "cube")
  (check-equal? (interface:universe.symmetry-transformations trunc-oct)
                (list (cons (array #[1/2 1/2 1/2])
                            (array #[#[1 0 0] #[0 1 0] #[ 0 0 1]]))))
  (check-equal? (interface:universe.convention trunc-oct) "my-own")
  (check-equal? (interface:universe.molecules trunc-oct)
                (list (cons (fragment "water" water) 10)))

  ; Generic utility functions

  (let ([da (fragment "da" di-alanine)])
    (check-equal? (sequence->list (interface:in-atoms da))
                  (append (interface:fragment.atoms
                           (interface:fragment.lookup-node da "ALA1"))
                          (interface:fragment.atoms
                           (interface:fragment.lookup-node da "ALA2"))))
    (check-equal? (sequence->list (interface:in-atoms-with-indices da))
                  (map list (sequence->list (interface:in-atoms da))
                            (sequence->list
                             (in-range (interface:number-of-atoms da)))))
    (check-equal? (sequence->list (interface:in-sites-with-indices da))
                  (map list (sequence->list (interface:in-atoms da))
                            (sequence->list
                             (in-range (interface:number-of-atoms da)))
                            (sequence->list
                             (in-range (interface:number-of-atoms da)))))
    (check-equal? (interface:number-of-atoms da) 23)
    (check-equal? (interface:number-of-sites da) 23))

  ; Dictionary interface

  (let ([w (fragment "w" water)])
    (check-equal? (dict-ref w "H1") (atom "H1" H 1))
    (check-equal? (dict-ref w "H2") (atom "H2" H 1))
    (check-equal? (dict-ref w "O") (atom "O" O 1)))

  (let ([da (fragment "da" di-alanine)])
    (check-equal? (dict-ref da "ALA1") (fragment "ALA1" alanine-nter))
    (check-equal? (dict-ref da "ALA2") (fragment "ALA2" alanine-cter))
    (check-equal? (dict-ref (dict-ref da "ALA1") "CA") (atom "CA" C 1)))

  ; Stream interface

  (let ([da (fragment "da" di-alanine)])
    (check-equal? (stream-first da) (fragment "ALA1" alanine-nter))
    (check-equal? (stream-first (stream-rest da))
                  (fragment "ALA2" alanine-cter))
    (check-equal? (stream-length da) 2))

  (let ([w (fragment "w" water)])
    (check-exn exn:fail:contract? (λ () (stream-length w)))
    (check-exn exn:fail:contract? (λ () (stream-first w)))
    (check-exn exn:fail:contract? (λ () (stream-rest w))))

  ; Factory function

  (check-equal? water10 (make-universe water10))
  (check-true (interface:universe.equivalent? water10 (make-universe water10))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Configurations
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(test-case "configuration"

  ; MOSAIC interface

  (let ([c (configuration water10 (make-array #(30 3) 0.)
                          (make-array #(0) (void)) 'float64)])
    (check-equal? (interface:configuration.universe c) water10)
    (check-equal? (interface:configuration.positions c)
                  (make-array #(30 3) 0.))
    (check-equal? (interface:configuration.cell-parameters c)
                  (make-array #(0) (void)))
    (check-equal? (interface:configuration.number-type c) "float64")))

; End of test-suite
))
