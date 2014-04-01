#lang racket

(require "interface.rkt"
         (prefix-in model: "model.rkt")
         math/array)

(provide messages
         validate-items
         validate-label
         validate-atom valid-atom-types valid-chemical-elements
         validate-fragment valid-polymer-types
         valid-bond-orders
         validate-universe
         valid-cell-shapes
         validate-configuration
         valid-cell-parameter-shapes
         valid-number-types-for-configurations)

; A bit of infrastructure code that permits writing the
; validation suites in a nice and compact way.

(define (mosaic-error text)
  (λ (path)
     (list (cons path text))))

(define no-error
  (void))

(define (combine . tests)
  (filter (negate void?) (flatten tests)))

(define (messages tests path)
  (apply append
         (map (λ (t) (t path))
                     (combine tests))))

(define (validate-node node . tests)
  (λ (path)
     (messages tests (cons node path))))

(define (validate-items items)
  (λ (path)
     (messages
        (for/list ([item items])
          (cond
           [(universe? item) (validate-universe item)]
           [else (validate-node item (mosaic-error "unknown item type"))]))
        path)))

; And now... the validation code.

; Labels are strings with certain conditions on their contents.
; They are used in various types of nodes.

(define valid-characters-in-labels 
  (list->set
   (string->list (string-append 
                  "abcdefghijklmnopqrstuvwxyz"
                  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                  "0123456789!#$%&?@^_~+-*/=,()[]'"))))

(define (validate-label-string l text)
  (cond
   [(not (string? l)) (mosaic-error (format "~a not a string: ~s" text l))]
   [(string=? l "") (mosaic-error (format "~a is empty" text))]
   [(> (string-length l) 32767) (mosaic-error (format "~a too long: ~a chars"
                                                      text (string-length l)))]
   [(let ([invalid-chars
           (filter (lambda (ch)
                     (not (set-member? valid-characters-in-labels ch)))
                   (string->list l))])
      (if (null? invalid-chars)
          no-error
          (mosaic-error (format "~a contains invalid character(s): ~s"
                                text (list->string invalid-chars)))))]
   [else no-error]))

(define (validate-label l)
  (validate-label-string l "label"))

; Atoms

(define (validate-atom a)
  (validate-node a
   (if (not (atom? a))
       (mosaic-error "not an atom")
       (combine
        (validate-label (node.label a))
        (validate-atom-type (atom.type a))
        (when (equal? (atom.type a) "element")
           (validate-chemical-element (atom.name a)))
        (validate-atom-nsites (atom.nsites a))))))

(define valid-atom-types
  (set "element" "cgparticle" "dummy" ""))

(define valid-chemical-elements
  (set "Ac" "Ag" "Al" "Am" "Ar" "As" "At" "Au"
       "B" "Ba" "Be" "Bh" "Bi" "Bk" "Br"
       "C" "Ca" "Cd" "Ce" "Cf" "Cl" "Cm"
       "Co" "Cn" "Cr" "Cs" "Cu"
       "D" "Db" "Ds" "Dy"
       "Er" "Es" "Eu"
       "F" "Fe" "Fm" "Fr"
       "Ga" "Gd" "Ge"
       "H" "He" "Hf" "Hg" "Ho" "Hs"
       "I" "In" "Ir"
       "K" "Kr"
       "La" "Li" "Lr" "Lu"
       "Md" "Mg" "Mn" "Mo" "Mt"
       "N" "Na" "Nb" "Nd" "Ne" "Ni" "No" "Np"
       "O" "Os"
       "P" "Pa" "Pb" "Pd" "Pm" "Po" "Pr" "Pt" "Pu"
       "Ra" "Rb" "Re" "Rf" "Rg" "Rh" "Rn" "Ru"
       "S" "Sb" "Sc" "Se" "Sg" "Si" "Sm" "Sn" "Sr"
       "Ta" "Tb" "Tc" "Te" "Th" "Ti" "Tl" "Tm"
       "U" "V" "W" "Xe" "Y" "Yb" "Zn" "Zr"))

(define (validate-atom-type t)
  (unless (set-member? valid-atom-types t)
    (mosaic-error (format "not a valid atom type: ~s" t))))

(define (validate-chemical-element n)
  (unless (set-member? valid-chemical-elements n)
    (mosaic-error (format "not a valid chemical element symbol: ~s" n))))

(define (validate-atom-nsites n)
  (unless (exact-positive-integer? n)
    (mosaic-error (format "not a positive integer: ~a" n))))


; Fragments

(define (validate-fragment f)
  (validate-node f
   (if (not (fragment? f))
       (mosaic-error "not a fragment")
       (combine
        (validate-label (node.label f))
        (validate-species (fragment.species f))
        (when (fragment.polymer? f)
            (validate-polymer-type (fragment.polymer-type f)))
        (for/list ([sf (fragment.subfragments f)])
          (validate-fragment sf))
        (for/list ([a (fragment.atoms f)])
          (validate-atom a))
        (for/list ([b (fragment.bonds f)])
          (validate-bond f b))
        (validate-labels-unique f)))))

(define valid-polymer-types
  (set ""
       "polypeptide"
       "polyribonucleotide"
       "polydeoxyribonucleotide"
       "polynucleotide"))

(define (validate-species l)
  (validate-label-string l "species name"))

(define (validate-polymer-type t)
  (unless (set-member? valid-polymer-types t)
    (mosaic-error (format "not a valid polymer type: ~a" t))))

(define valid-bond-orders
  (set ""
       "single"
       "double"
       "triple"
       "quadruple"
       "aromatic"))

(define (validate-bond fragment bond)
  (match bond
    [(vector p1 p2 order)
     (combine
      (validate-path fragment p1)
      (validate-path fragment p2)
      (unless (set-member? valid-bond-orders order)
        (mosaic-error (format "not a valid bond order: ~s" order))))]))

(define (validate-path fragment path)
  (let ([node-path (fragment-lookup-path fragment path)])
    (cond
     [(not node-path)
         (mosaic-error (format "no atom ~a" path))]
     [(not (atom? (first node-path)))
         (mosaic-error (format "node ~a is not an atom" path))])))

(define (validate-labels-unique f)
  (let* ([sfs (fragment.subfragments f)]
         [sf-labels (apply set (map node.label sfs))]
         [as  (fragment.atoms f)]
         [a-labels (apply set (map node.label as))])
    (unless (= (set-count (set-union sf-labels a-labels))
                    (+ (length sfs) (length as)))
           (mosaic-error "duplicate label(s)"))))


; Universes

(define (validate-universe u)
  (validate-node u
   (if (not (universe? u))
       (mosaic-error "not a universe")
       (combine
        (validate-cell-shape (universe.cell-shape u))
        (validate-symmetry-transformations
           (universe.symmetry-transformations u))
        (validate-convention (universe.convention u))
        (for/list ([m (universe.molecules u)])
          (match m
            [(cons fragment count)
               (combine (validate-fragment fragment)
                 (unless (exact-positive-integer? count)
                   (mosaic-error (format "not a positive integer: ~s" count))))]
            [term
               (mosaic-error
                (format "not a fragment-count pair: ~a" term))]))))))

(define valid-cell-shapes
  (set "infinite"
       "cube"
       "cuboid"
       "parallelepiped"))

(define (validate-cell-shape cs)
  (unless (set-member? valid-cell-shapes cs)
    (mosaic-error (format "not a valid cell shape: ~s" cs))))

(define (validate-convention c)
  (validate-label-string c "convention"))

(define (validate-symmetry-transformations st)
  (combine
   (for/list ([t st])
     (if (pair? t)
         (match t
           [(cons trans rot)
            (combine
             (validate-translation trans)
             (validate-rotation rot))])
         (mosaic-error (format "not a pair: ~a" t))))))

(define (validate-translation t)
  (unless (and (array? t)
               (equal? (array-shape t) #[3]))
    (mosaic-error (format "not an array of shape #[3]: ~a" t))))

(define (validate-rotation r)
  (unless (and (array? r)
               (equal? (array-shape r) #[3 3]))
    (mosaic-error (format "not an array of shape #[3 3]: ~a" r))))

; Configurations

(define (validate-configuration c)
  (validate-node c
    (if (not (configuration? c))
        (mosaic-error "not a configuration")
        (combine
         (validate-number-type-conf (configuration.number-type c))
         (let ([u (configuration.universe c)])
           (if (not (universe? u))
               (mosaic-error (format "not a universe: ~a" u))
               (combine
                (validate-universe u)
                (validate-positions (configuration.positions c)
                                    (number-of-sites u))
                (validate-cell-parameters (configuration.cell-parameters c)
                                          (universe.cell-shape u)))))))))

(define valid-cell-parameter-shapes
  #hash(("infinite" . #[0])
        ("cube" . #[])
        ("cuboid" . #[3])
        ("parallelepiped" . #[3 3])))

(define (validate-cell-parameters cp cs)
  (if (not (array? cp))
      (mosaic-error (format "not an array: ~a" cp))
      (let ([cp-shape (array-shape cp)]
            [ucp-shape (dict-ref valid-cell-parameter-shapes cs)])
        (unless (equal? cp-shape ucp-shape)
          (mosaic-error (format "invalid cell parameter shape ~a (should be ~a)"
                                cp-shape ucp-shape ))))))

(define valid-number-types-for-configurations
  (set "float32" "float64"))

(define (validate-number-type-conf nt)
  (unless (set-member? valid-number-types-for-configurations nt)
    (mosaic-error (format "invalid number type: ~a" nt))))

(define (validate-positions pos nsites)
  (if (not (array? pos))
      (mosaic-error (format "not an array: ~a" pos))
      (let ([pos-shape (array-shape pos)]
            [upos-shape (vector nsites 3)])
        (unless (equal? pos-shape upos-shape)
          (mosaic-error (format "invalid position shape ~a (should be ~a)"
                                pos-shape upos-shape))))))
