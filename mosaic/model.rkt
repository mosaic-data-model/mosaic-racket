#lang racket

(provide (all-defined-out))

(require (prefix-in interface:  "interface.rkt")
         (for-syntax syntax/parse)
         math/array
         generic-bind)

; Atom and Fragment descriptors

(struct atom-descr (type name)
    #:transparent)

(define (~frag-wrapper fn)
  (lambda (fragment-descr . args)
    (apply fn (fragment "_" fragment-descr) args)))

(struct fragment-descr (species subfragments atoms bonds)
    #:transparent
    #:methods gen:dict
    [(define dict-ref
       (~frag-wrapper interface:fragment-dict-ref))
     (define dict-count
       (~frag-wrapper interface:fragment-dict-count))
     (define dict-iterate-first
       (~frag-wrapper interface:fragment-dict-iterate-first))
     (define dict-iterate-next
       (~frag-wrapper interface:fragment-dict-iterate-next))
     (define dict-iterate-key
       (~frag-wrapper interface:fragment-dict-iterate-key))
     (define dict-iterate-value
       (~frag-wrapper interface:fragment-dict-iterate-value))])

(struct polymer-descr (type species subfragments bonds)
    #:transparent
    #:methods gen:dict
    [(define dict-ref
       (~frag-wrapper interface:fragment-dict-ref))
     (define dict-count
       (~frag-wrapper interface:fragment-dict-count))
     (define dict-iterate-first
       (~frag-wrapper interface:fragment-dict-iterate-first))
     (define dict-iterate-next
       (~frag-wrapper interface:fragment-dict-iterate-next))
     (define dict-iterate-key
       (~frag-wrapper interface:fragment-dict-iterate-key))
     (define dict-iterate-value
       (~frag-wrapper interface:fragment-dict-iterate-value))])

; Atoms and fragments

(struct node (label)
    #:transparent
    #:methods interface:gen:node
    [(define (node.label atom)
       (node-label atom))])

(struct atom node (adescr nsites)
    #:transparent
    #:methods interface:gen:atom
    [(define (atom.type atom)
       (let ([t (atom-descr-type (atom-adescr atom))])
         (if (equal? t 'unknown)
             ""
             (symbol->string t))))
     (define (atom.name atom)
       (atom-descr-name (atom-adescr atom)))
     (define (atom.nsites atom)
       (atom-nsites atom))]
    #:methods gen:custom-write
    [(define (write-proc atom port mode)
       (write-string "<atom " port)
       (write-string (node-label atom) port)
       (write-string ">" port))])

(define (make-atom a)
  (atom (interface:node.label a)
        (atom-descr (let ([t (interface:atom.type a)])
                      (if (string=? t "") 'unknown (string->symbol t))) 
                    (interface:atom.name a))
        (interface:atom.nsites a)))

(struct fragment node (fdescr)
    #:transparent
    #:methods interface:gen:fragment
    [(define (fragment.species fragment)
       (let ([fd (fragment-fdescr fragment)])
         (if (polymer-descr? fd)
             (polymer-descr-species fd)
             (fragment-descr-species fd))))
     (define (fragment.subfragments fragment)
       (let ([fd (fragment-fdescr fragment)])
         (if (polymer-descr? fd)
             (polymer-descr-subfragments fd)
             (fragment-descr-subfragments fd))))
     (define (fragment.atoms fragment)
       (let ([fd (fragment-fdescr fragment)])
         (if (polymer-descr? fd)
             '()
             (fragment-descr-atoms fd))))
     (define fragment.lookup-node interface:fragment-lookup-node)
     (define (fragment.bonds fragment)
       (map (lambda (b)
              (match-let* ([(list a1 a2) (set->list (car b))]
                           [order (cdr b)]
                           [m-order (if (equal? order 'unknown)
                                        ""
                                        (symbol->string order))])
                (if (string<? a1 a2)
                    (vector-immutable a1 a2 m-order)
                    (vector-immutable a2 a1 m-order))))
            (set->list (let ([fd (fragment-fdescr fragment)])
                         (if (polymer-descr? fd)
                             (polymer-descr-bonds fd)
                             (fragment-descr-bonds fd))))))
     (define (fragment.polymer? fragment)
       (polymer-descr? (fragment-fdescr fragment)))
     (define (fragment.polymer-type fragment)
       (symbol->string (polymer-descr-type (~get-polymer-descr fragment))))]
    #:methods gen:dict
    [(define dict-ref interface:fragment-dict-ref)
     (define dict-count interface:fragment-dict-count)
     (define dict-iterate-first interface:fragment-dict-iterate-first)
     (define dict-iterate-next interface:fragment-dict-iterate-next)
     (define dict-iterate-key interface:fragment-dict-iterate-key)
     (define dict-iterate-value interface:fragment-dict-iterate-value)]
    #:methods gen:stream
    [(define (stream-empty? fragment)
       (empty? (polymer-descr-subfragments (~get-polymer-descr fragment))))
     (define (stream-first fragment)
       (first (polymer-descr-subfragments (~get-polymer-descr fragment))))
     (define (stream-rest fragment)
       (rest (polymer-descr-subfragments (~get-polymer-descr fragment))))]
    #:methods gen:custom-write
    [(define (write-proc fragment port mode)
       (write-string (if (polymer-descr? (fragment-fdescr fragment))
                         "<polymer "
                         "<fragment ") port)
       (write-string (node-label fragment) port)
       (write-string ">" port))])

(define (~get-polymer-descr fragment)
  (let ([fd (fragment-fdescr fragment)])
    (if (polymer-descr? fd)
        fd
        (raise-argument-error 'fragment.polymer-type
                              "fragment.polymer?"
                              fragment))))

(define (make-fragment f)
  (let ([species (interface:fragment.species f)]
        [subfragments (map make-fragment
                           (interface:fragment.subfragments f))]
        [atoms (map make-atom (interface:fragment.atoms f))]
        [bonds (list->set
                (map (λ (b)
                        (cons (set (vector-ref b 0)
                                   (vector-ref b 1))
                              (string->symbol (vector-ref b 2))))
                     (interface:fragment.bonds f)))])
    (fragment (interface:node.label f)
              (if (interface:fragment.polymer? f)
                  (polymer-descr (string->symbol
                                  (interface:fragment.polymer-type f))
                                 species subfragments bonds)
                  (fragment-descr species subfragments atoms bonds)))))

; Universe

(struct universe (cell-shape symmetry-transformations convention molecules)
    #:transparent
    #:methods interface:gen:universe
    [(define (universe.cell-shape universe)
       (symbol->string (universe-cell-shape universe)))
     (define (universe.symmetry-transformations universe)
       (universe-symmetry-transformations universe))
     (define (universe.convention universe)
       (universe-convention universe))
     (define (universe.molecules universe)
       (universe-molecules universe))]
    #:methods gen:custom-write
    [(define (write-proc universe port mode)
       (write-string "<" port)
       (write-string (symbol->string (universe-cell-shape universe)) port)
       (write-string " universe>" port))])

(define (make-universe u)
  (universe (string->symbol (interface:universe.cell-shape u))
            (interface:universe.symmetry-transformations u)
            (interface:universe.convention u)
            (~for/list ([($: fragment count)
                         (interface:universe.molecules u)])
               (cons (make-fragment fragment) count))))

; Macros for conveniently defining atoms and fragments

(define-syntax-rule (define-element symbol)
  (define symbol (atom-descr 'element (symbol->string (quote symbol)))))

(define-syntax define-elements
  (syntax-rules ()
    [(define-elements s)
     (define s (atom-descr 'element (symbol->string (quote s))))]
    [(define-elements s1 s2 ...)
     (begin (define-element s1)
            (define-elements s2 ...))]))

(define-elements 
  Ac Ag Al Am Ar As At Au
  B Ba Be Bh Bi Bk Br
  C Ca Cd Ce Cf Cl Cm
  Co Cn Cr Cs Cu
  D Db Ds Dy
  Er Es Eu
  F Fe Fm Fr
  Ga Gd Ge
  H He Hf Hg Ho Hs
  I In Ir
  K Kr
  La Li Lr Lu
  Md Mg Mn Mo Mt
  N Na Nb Nd Ne Ni No Np
  O Os
  P Pa Pb Pd Pm Po Pr Pt Pu
  Ra Rb Re Rf Rg Rh Rn Ru
  S Sb Sc Se Sg Si Sm Sn Sr
  Ta Tb Tc Te Th Ti Tl Tm
  U V W Xe Y Yb Zn Zr)


(begin-for-syntax

  (define-syntax-class atom-spec
   #:description "atom specification"
   (pattern (label:expr atom-descr:expr)
            #:with spec #'(atom label atom-descr 1))
   (pattern (label:expr atom-descr:expr nsites:expr)
            #:with spec #'(atom label atom-descr nsites)))

 (define-syntax-class fragment-spec
   #:description "fragment specification"
   (pattern (label:expr fragment-descr:expr)
            #:with spec #'(fragment label fragment-descr)))

 (define-syntax-class bond-spec
   #:description "bond specification"
   (pattern (atom1:expr atom2:expr)
            #:with spec #'(cons (set atom1 atom2) 'unknown))
   (pattern (atom1:expr atom2:expr order:expr)
            #:with spec #'(cons (set atom1 atom2) order))))

(define-syntax (define-fragment stx)
  (syntax-parse stx
   [(define-fragment species
          (~optional (~seq #:subfragments (frag:fragment-spec ...)))
          (~optional (~seq #:atoms (atom:atom-spec ...)))
          (~optional (~seq #:bonds (bond:bond-spec ...))))
    (with-syntax ([species-str  (symbol->string (syntax->datum #'species))]
                  [frag (or (attribute frag.spec) '())]
                  [atom (or (attribute atom.spec) '())]
                  [bond (or (attribute bond.spec) '())])
      #'(define species
          (fragment-descr species-str
                          (list . frag) (list . atom) (set . bond))))]))

(define-syntax (define-polymer stx)
  (syntax-parse stx
   [(define-polymer species polymer-type
          (~optional (~seq #:subfragments (frag:fragment-spec ...)))
          (~optional (~seq #:bonds (bond:bond-spec ...))))
    (with-syntax ([species-str  (symbol->string (syntax->datum #'species))]
                  [frag (or (attribute frag.spec) '())]
                  [bond (or (attribute bond.spec) '())])
      #'(define species
          (polymer-descr polymer-type species-str
                         (list . frag) (set . bond))))]))

; Configurations

(struct configuration (universe positions cell-parameters number-type)
    #:transparent
    #:methods interface:gen:configuration
    [(define (configuration.universe configuration)
       (configuration-universe configuration))
     (define (configuration.positions configuration)
       (configuration-positions configuration))
     (define (configuration.cell-parameters configuration)
       (configuration-cell-parameters configuration))
     (define (configuration.number-type configuration)
       (symbol->string (configuration-number-type configuration)))])

(define (make-configuration c)
  (configuration (make-universe (interface:configuration.universe c))
                 (interface:configuration.positions c)
                 (interface:configuration.cell-parameters c)
                 (string->symbol (interface:configuration.number-type c))))

; Generic factory function

(define (make-data-item d)
  (cond
   [(interface:universe? d) (make-universe d)]
   [(interface:configuration? d) (make-configuration d)]
   [else (error "not a Mosaic data item")]))
