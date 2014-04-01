#lang racket

(provide (all-defined-out))

(require racket/generic
         racket/stream
         racket/generator
         generic-bind
         (for-syntax syntax/parse))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Version information
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define mosaic-version-major 1)
(define mosaic-version-minor 0)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Support code for comparisons
;
; Comparisons are implemented as generators that
; return at each step one difference between the
; two data structures. These generators are the
; basis for two functions that implement the two
; most common use cases: checking for equivalence
; (no differences), and returning a complete list
; of differences.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-for-syntax (append-id stx x y)
  (define x:str (if (string? x) x (symbol->string (syntax-e x))))
  (define y:str (if (string? y) y (symbol->string (syntax-e y))))
  (datum->syntax stx (string->symbol (string-append x:str y:str))))

(define-syntax (define-comparison stx)
  (syntax-parse stx
    [(_ kind:id fn:expr)
     (define compare-name (append-id stx #'kind ".compare"))
     (define equivalent?-name (append-id stx #'kind ".equivalent?"))
     (define diffs-name (append-id stx #'kind ".diffs"))
     #`(begin
         (define #,compare-name fn)
         (define (#,equivalent?-name v1 v2)
           (stream-empty? (sequence->stream (#,compare-name v1 v2))))
         (define (#,diffs-name v1 v2)
           (sequence->list (#,compare-name v1 v2))))]))

(define-syntax-rule (cmp label predicate v1 v2)
  (cond [(~cmp label predicate v1 v2) => yield]))

(define (~cmp label predicate v1 v2)
  (if (predicate v1 v2)
      #f
      (list label v1 v2)))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Universe interface
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-generics universe
  [universe.cell-shape universe]
  [universe.symmetry-transformations universe]
  [universe.convention universe]
  [universe.molecules universe])

(define-generics node
  [node.label node])

(define-generics fragment
  [fragment.species fragment]
  [fragment.subfragments fragment]
  [fragment.atoms fragment]
  [fragment.lookup-node fragment label]
  [fragment.bonds fragment]
  [fragment.polymer? fragment]
  [fragment.polymer-type fragment])

(define-generics atom
  [atom.type atom]
  [atom.name atom]
  [atom.nsites atom])

(define-comparison universe
  (λ  (u1 u2)
      (in-generator
       (cmp 'universe.cell-shape
            equal? (universe.cell-shape u1) (universe.cell-shape u2))
       (cmp 'universe.convention
            equal? (universe.convention u1) (universe.convention u2))
       (cmp 'universe.symmetry-transformations
            equal? (sequence->list (universe.symmetry-transformations u1))
            (sequence->list (universe.symmetry-transformations u2)))
       (let ([ms1 (sequence->list (universe.molecules u1))]
             [ms2 (sequence->list (universe.molecules u2))])
         (cmp 'universe.number-of-molecules = (length ms1) (length ms2))
         (~for ([($: f1 c1) ms1]
                [($: f2 c2) ms2])
               (cmp 'universe.molecule-count = c1 c2)
               (for ([fcmp (fragment.compare f1 f2)]) (yield fcmp)))))))

(define-comparison fragment
  (λ (f1 f2)
     (in-generator
      (cmp 'node.label equal? (node.label f1) (node.label f2))
      (cmp 'fragment.species equal?
           (fragment.species f1) (fragment.species f2))
      (cmp 'fragment.polymer? equal?
           (fragment.polymer? f1) (fragment.polymer? f2))
      (when (fragment.polymer? f1)
        (cmp 'fragment.polymer-type equal?
             (fragment.polymer-type f1) (fragment.polymer-type f2)))
      (let ([fs1 (sequence->list (fragment.subfragments f1))]
            [fs2 (sequence->list (fragment.subfragments f2))])
        (cmp 'fragment.number-of-subfragments = (length fs1) (length fs2))
        (for ([f1 fs1]
              [f2 fs2])
          (for ([fcmp (fragment.compare f1 f2)]) (yield fcmp))))
      (let ([as1 (sequence->list (fragment.atoms f1))]
            [as2 (sequence->list (fragment.atoms f2))])
        (cmp 'fragment.number-of-atoms = (length as1) (length as2))
        (for ([a1 as1]
              [a2 as2])
          (for ([acmp (atom.compare a1 a2)]) (yield acmp))))
      (let ([bs1 (for/set ([b (fragment.bonds f1)])
                   (match b [(vector a1 a2 order) (cons (set a1 a2) order)]))]
            [bs2 (for/set ([b (fragment.bonds f2)])
                   (match b [(vector a1 a2 order) (cons (set a1 a2) order)]))])
        (cmp 'fragment.number-of-bonds = (set-count bs1) (set-count bs2))
        (when (= (set-count bs1) (set-count bs2))
          (cmp 'fragment.bonds equal? bs1 bs2))))))

(define-comparison atom
  (λ (a1 a2)
     (in-generator
      (cmp 'node.label equal? (node.label a1) (node.label a2))
      (cmp 'atom.type equal? (atom.type a1) (atom.type a2))
      (cmp 'atom.name equal? (atom.name a1) (atom.name a2))
      (cmp 'atom.nsites = (atom.nsites a1) (atom.nsites a2)))))


; Node lookup by label and path.

; fragment-lookup-node is provided as a default implementation
; for fragment.lookup-node.
(define (fragment-lookup-node fragment label)
  (or (for/first ([sf (fragment.subfragments fragment)]
                  #:when (equal? (node.label sf) label))
        sf)
      (for/first ([a (fragment.atoms fragment)]
                  #:when (equal? (node.label a) label))
        a)))

; fragment-lookup-path performs sequential lookup for
; the elements in a path, and returns the reverse of the list
; of traversed nodes.
(define (fragment-lookup-path fragment path)
  (let loop ([node-path (list fragment)]
             [split-path (string-split path "." #:trim? #f)])
    (let ([node (first node-path)])
      (cond
       [(not node)  #f]
       [(empty? split-path)  node-path]
       [(not (fragment? node))  #f]
       [else (loop (cons (fragment.lookup-node node (first split-path))
                         node-path)
                   (rest split-path))]))))

; A default implementation for Racket's dict interface for fragments,
; based only on the interface functions. For every lookup, it iterates
; over the subnodes. This is not as bad as it sounds because the number
; of subnodes is small in most fragment definitions.

(struct ~dict-index (pos sf-count atom-count))

(define (fragment-dict-ref fragment key
                  [default (lambda () (error "key not found" key))])
  (or (and (string? key)
           (fragment.lookup-node fragment key))
      (if (procedure? default) (default) default)))

(define (fragment-dict-count fragment)
  (+ (length (fragment.subfragments fragment))
     (length (fragment.atoms fragment))))

(define (fragment-dict-iterate-first fragment)
  (if (and (empty? (fragment.subfragments fragment))
           (empty? (fragment.atoms fragment))) 
      #f
      (~dict-index 0
                   (length (fragment.subfragments fragment))
                   (length (fragment.atoms fragment)))))

(define (fragment-dict-iterate-next fragment index)
  (when (not (~dict-index? index))
    (raise (exn:fail:contract "invalid index"
                              (current-continuation-marks))))
  (let ([next (add1 (~dict-index-pos index))])
    (if (equal? next (+ (~dict-index-sf-count index)
                        (~dict-index-atom-count index)))
        #f
        (~dict-index next
                     (~dict-index-sf-count index)
                     (~dict-index-atom-count index)))))

(define (fragment-dict-iterate-value fragment index)
  (when (not (~dict-index? index))
    (raise (exn:fail:contract "invalid index"
                              (current-continuation-marks))))
  (let* ([nsf (~dict-index-sf-count index)]
         [pos (~dict-index-pos index)])
    (if (< pos nsf)
        (list-ref (fragment.subfragments fragment) pos)
        (list-ref (fragment.atoms fragment) (- pos nsf)))))

(define (fragment-dict-iterate-key fragment pos)
  (node.label (fragment-dict-iterate-value fragment pos)))

; Utility functions defined in terms of the interface functions

(define (atom-sequence ms)
  (define (fragment-atom-sequence fragment)
    (stream-append
      (apply stream-append
             (map fragment-atom-sequence
                  (fragment.subfragments fragment)))
      (fragment.atoms fragment)))
  (define (universe-fragment-sequence universe)
    (in-generator
     (~for ([($: f c) (universe.molecules universe)])
        (for ([_ (in-range c)])
          (yield f)))))
  (define (universe-atom-sequence universe)
    (apply stream-append
      (map fragment-atom-sequence
           (sequence->list (universe-fragment-sequence universe)))))
  (cond
   [(atom? ms) (list ms)]
   [(fragment? ms) (fragment-atom-sequence ms)]
   [(universe? ms) (universe-atom-sequence ms)]))

(define (number-of-atoms ms)
  (sequence-length (atom-sequence ms)))

(define (number-of-sites ms)
  (sequence-fold + 0 (sequence-map atom.nsites (atom-sequence ms))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Configuration interface
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-generics configuration
  [configuration.universe configuration]
  [configuration.positions configuration]
  [configuration.cell-parameters configuration]
  [configuration.number-type configuration])

(define-comparison configuration
  (λ (c1 c2)
     (in-generator
      (cmp 'configuration.universe
           universe.equivalent? (configuration.universe c1)
           (configuration.universe c2))
      (cmp 'configuration.positions
           equal? (configuration.positions c1) (configuration.positions c2))
      (cmp 'configuration.cell-parameters
           equal? (configuration.cell-parameters c1)
           (configuration.cell-parameters c2))
      (cmp 'configuration.number-type equal?
           (configuration.number-type c1) (configuration.number-type c2)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Generic functions for any data item type
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (compare item1 item2)
  (cond
   [(and (universe? item1)
         (universe? item2)) (universe.compare item1 item2)]
   [(and (configuration? item1)
         (configuration? item2)) (configuration.compare item1 item2)]
   [else (list 'incompatible-items item1 item2)]))

(define (equivalent? v1 v2)
  (stream-empty? (sequence->stream (compare v1 v2))))

(define (diffs v1 v2)
  (sequence->list (compare v1 v2)))
