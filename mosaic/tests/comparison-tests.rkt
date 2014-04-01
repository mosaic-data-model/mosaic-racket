#lang racket

(require rackunit
         math/array
          "../interface.rkt"
         (prefix-in model: "../model.rkt"))

(provide comparison-tests)

; Fragment definition used in the tests
(model:define-fragment water
                       #:atoms
                       (["H1" model:H]
                        ["H2" model:H]
                        ["O" model:O])
                       #:bonds
                       (["H1" "O" 'single]
                        ["H2" "O" 'single]))

(define comparison-tests
  (test-suite "comparison"

; Atoms
(test-case "atoms"

  (check-true (atom.equivalent? (model:atom "C" model:C 1)
                                (model:atom "C" model:C 1)))

  (let ([a1 (model:atom "C1" model:C 1)]
        [a2 (model:atom "C2" model:C 1)])
    (check-equal? (atom.diffs a1 a2)
                  '((node.label "C1" "C2"))))

  (let ([a1 (model:atom "H" (model:atom-descr 'element "H") 1)]
        [a2 (model:atom "H" (model:atom-descr 'dummy "H") 1)])
    (check-equal? (atom.diffs a1 a2)
                  '((atom.type "element" "dummy"))))

  (let ([a1 (model:atom "H" (model:atom-descr 'element "H") 1)]
        [a2 (model:atom "H" (model:atom-descr 'element "O") 1)])
    (check-equal? (atom.diffs a1 a2)
                  '((atom.name "H" "O"))))

  (let ([a1 (model:atom "C" model:C 1)]
        [a2 (model:atom "C" model:C 2)])
    (check-equal? (atom.diffs a1 a2)
                  '((atom.nsites 1 2)))))

; Fragments
(test-case "fragments"

  (model:define-fragment water-permuted-bonds
                         #:atoms
                         (["H1" model:H]
                          ["H2" model:H]
                          ["O" model:O])
                         #:bonds
                         (["O" "H2" 'single]
                          ["H1" "O" 'single]))

  (model:define-fragment water-three-bonds
                         #:atoms
                         (["H_one" model:H]
                          ["H_two" model:H]
                          ["O" model:O])
                         #:bonds
                         (["H_one" "O" 'single]
                          ["H_two" "O" 'single]
                          ["H_one" "H_two" 'single]))

  (check-true (fragment.equivalent? (model:fragment "w" water)
                                    (model:fragment "w" water)))

  (let ([f1 (model:fragment "w1" water)]
        [f2 (model:fragment "w2" water)])
    (check-equal? (fragment.diffs f1 f2)
                  '((node.label "w1" "w2"))))

  (let ([f1 (model:fragment "w" water)]
        [f2 (model:fragment "w" water-three-bonds)])
    (check-equal? (fragment.diffs f1 f2)
                  '((fragment.species "water" "water-three-bonds")
                    (node.label "H1" "H_one")
                    (node.label "H2" "H_two")
                    (fragment.number-of-bonds 2 3))))

  (let ([f1 (model:fragment "w" water)]
        [f2 (model:fragment "w" water-permuted-bonds)])
    (check-equal? (fragment.diffs f1 f2)
                  '((fragment.species "water" "water-permuted-bonds")))))

; Universes
(test-case "universe"

  (define water-infinite
    (model:universe 'infinite '() "my-own"
                    (list (cons (model:fragment "water" water)  10))))

  (define water-trunc-oct
    (model:universe 'cube
                    (list (cons (array #[1/2 1/2 1/2])
                                (array #[#[1 0 0] #[0 1 0] #[ 0 0 1]])))
                    "my-own"
                    (list (cons (model:fragment "water" water) 10))))

  (check-true (universe.equivalent? water-infinite water-infinite))
  (check-true (universe.equivalent? water-trunc-oct water-trunc-oct))

  (check-equal? (universe.diffs water-infinite water-trunc-oct)
                (list '(universe.cell-shape "infinite" "cube")
                      (list 'universe.symmetry-transformations
                            '()
                            (list (cons (array #[1/2 1/2 1/2])
                                        (array #[#[1 0 0] #[0 1 0] #[ 0 0 1]])))))))))
