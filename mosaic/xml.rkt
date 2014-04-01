#lang racket

(provide items-from-xml sxml->mosaic sxml-items
         mosaic->sxml mosaic-sequence->sxml-sequence items-to-xml)

(require sxml
         math/array
         racket/generator
         (prefix-in interface: "interface.rkt")
         generic-bind)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read Mosaic data from XML
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; A stream interface for the data structure returned by the lazy SSAX
; routines, which are lists whose last element is a promise yielding
; the rest of the list.
(struct sxml-stream (results)
  #:methods gen:stream
  [(define (stream-empty? stream)
     (let ([r (sxml-stream-results stream)])
       (or (empty? r)
           (and (promise? (first r))
                (empty? (force (first r)))))))
   (define (stream-first stream)
     (let* ([r (sxml-stream-results stream)]
            [f (first r)])
       (lazy:node->sxml (if (promise? f)
                            (first (force f))
                            f))))
   (define (stream-rest stream)
     (let* ([r (sxml-stream-results stream)]
            [f (first r)])
       (sxml-stream (if (promise? f)
                        (rest (force f))
                        (rest r)))))])

; The base struct type for Mosaic SXML data
(struct mosaic-sxml (node) #:transparent)

; Return a stream of sxml nodes describing Mosaic items
(define (sxml-items port [xpath "//mosaic/*"])
  (sxml-stream ((lazy:sxpath xpath)
                (lazy:xml->sxml port '()))))

; Return a pair (xml-id . mosaic-item) for an SXML node
(define (sxml->mosaic node [refs #hash()])
  (let ([xml-id (sxml:attr node 'id)])
    (cons xml-id
          (match (sxml:name node)
            ['universe (sxml-universe node)]
            ['configuration (sxml-configuration node refs)]
            [else (sxml-mosaic-item node)]))))

; Return a sequence of (xml-id . mosaic-item)
(define (items-from-xml port [xpath "//mosaic/*"])
  (in-generator
   (let ([refs (make-hash)])
     (for ([node (sxml-items port xpath)])
       (match-let ([(cons xml-id item) (sxml->mosaic node refs)])
         (hash-set! refs xml-id item)
         (yield (cons xml-id item)))))))

; Extract data from unique subnodes
(define (subelement tag node)
  (let ([elements (filter (λ (n) (equal? (first n) tag))
                          (sxml:content node))])
    (unless (= 1 (length elements))
      (error (format "element ~a must contain exactly one ~a"
                     (first node) tag)))
    (first elements)))

(define (text-of tag node)
  (sxml:text (subelement tag node)))

(define (text->vector text)
  (list->vector (map string->number (string-split text))))

(define (text->array text [element-shape #[]] #:flonum [fl-flag #f])
  (let* ([linear (list->array
                  (map string->number (string-split text)))]
         [linear (if fl-flag (array->flarray linear) linear)]
         [n (array-size linear)]
         [el (array-size (make-array element-shape (void)))]
         [n-el (/ n el)]
         [shape (vector-append (vector n-el) element-shape)])
    (unless (exact-integer? n-el)
      (error (format "~a not a multiple of ~a" n el)))
    (array-reshape linear shape)))

(define (array-of tag node [element-shape #[]] #:flonum [fl-flag #f])
  (text->array (text-of tag node) element-shape))

; Convert arrays to strings
(define (array->text a)
  (values
   (string-join (map number->string (vector->list (array-shape a))))
   (string-join (map number->string (array->list (array-flatten a))))))

; Mosaic interface for universe SXML nodes
(struct sxml-universe mosaic-sxml ()  #:transparent
    #:methods interface:gen:universe
    [(define (universe.cell-shape universe)
       (sxml:attr (mosaic-sxml-node universe) 'cell_shape))
     (define (universe.symmetry-transformations universe)
       (let ([st ((sxpath "symmetry_transformations/transformation")
                  (mosaic-sxml-node universe))])
         (map (λ (t)
                 (let ([trans (array-of 'translation t #:flonum #t)]
                       [rot (array-of 'rotation t #[3] #:flonum #t)])
                   (cons trans rot)))
              st)))
     (define (universe.convention universe)
       (sxml:attr (mosaic-sxml-node universe) 'convention))
     (define (universe.molecules universe)
       (map (λ (molecule)
               (cons (sxml-fragment (subelement 'fragment molecule))
                     (sxml:num-attr molecule 'count)))
            ((sxpath "molecules/molecule") (mosaic-sxml-node universe))))])

; Mosaic interface for fragment SXML nodes
(struct sxml-fragment mosaic-sxml () #:transparent
    #:methods interface:gen:node
    [(define (node.label fragment)
       (sxml:attr (mosaic-sxml-node fragment) 'label))]
    #:methods interface:gen:fragment
    [(define (fragment.species fragment)
       (sxml:attr (mosaic-sxml-node fragment) 'species))
     (define (fragment.subfragments fragment)
       (map sxml-fragment
            ((sxpath "fragments/fragment") (mosaic-sxml-node fragment))))
     (define (fragment.atoms fragment)
       (map sxml-atom
            ((sxpath "atoms/atom") (mosaic-sxml-node fragment))))
     (define fragment.lookup-node interface:fragment-lookup-node)
     (define (fragment.bonds fragment)
       (map (λ (bond)
               (match-let ([(list a1 a2) (string-split (sxml:attr bond 'atoms))])
                 (vector a1 a2 (sxml:attr bond 'order))))
            ((sxpath "bonds/bond") (mosaic-sxml-node fragment))))
     (define (fragment.polymer? fragment)
       (not (not (sxml:attr (mosaic-sxml-node fragment) 'polymer_type))))
     (define (fragment.polymer-type fragment)
       (sxml:attr (mosaic-sxml-node fragment) 'polymer_type))])

; Mosaic interface for atom SXML nodes
(struct sxml-atom mosaic-sxml () #:transparent
    #:methods interface:gen:node
    [(define (node.label atom)
       (sxml:attr (mosaic-sxml-node atom) 'label))]
    #:methods interface:gen:atom
    [(define (atom.type atom)
       (sxml:attr (mosaic-sxml-node atom) 'type))
     (define (atom.name atom)
       (sxml:attr (mosaic-sxml-node atom) 'name))
     (define (atom.nsites atom)
       (let ([ns (sxml:num-attr (mosaic-sxml-node atom) 'nsites)])
         (if ns ns 1)))])

; Mosaic interface for configuration SXML nodes
(struct sxml-configuration mosaic-sxml (refs) #:transparent
    #:methods interface:gen:configuration
    [(define (configuration.universe configuration)
       (let* ([u-node (subelement 'universe (mosaic-sxml-node configuration))]
              [u-ref (sxml:attr u-node 'ref)]
              [refs (sxml-configuration-refs configuration)])
         (if u-ref
             (dict-ref refs u-ref)
             (sxml-universe u-node))))
     (define (configuration.positions configuration)
       (array-of 'positions (mosaic-sxml-node configuration) #[3] #:flonum #t))
     (define (configuration.cell-parameters configuration)
       (let* ([cp (subelement 'cell_parameters
                              (mosaic-sxml-node configuration))]
              [cp-shape (text->vector (sxml:attr cp 'shape))])
         (array-reshape (text->array (sxml:text cp) #:flonum #t) cp-shape)))
     (define (configuration.number-type configuration)
       (sxml:attr (subelement 'positions (mosaic-sxml-node configuration))
                  'type))])

; Temporary catch-all for unimplemented item types
(struct sxml-mosaic-item mosaic-sxml () #:transparent)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Write Mosaic data to XML
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (sxml-from-universe universe xml-id refs)

  (define (sxml-from-molecule m)
    `(molecule (@ (count ,(number->string (cdr m))))
               ,(sxml-from-fragment (car m))))

  (define (sxml-from-fragment f)
    (let ([subfragments (map sxml-from-fragment
                             (interface:fragment.subfragments f))]
          [atoms (map sxml-from-atom
                      (interface:fragment.atoms f))]
          [bonds (map sxml-from-bond
                      (interface:fragment.bonds f))])
      `(fragment (@ (label ,(interface:node.label f))
                    (species ,(interface:fragment.species f))
                    ,@(if (interface:fragment.polymer? f)
                          `((polymer_type ,(interface:fragment.polymer-type f)))
                          '()))
                 ,@(if (empty? subfragments) '() `((fragments ,@subfragments)))
                 ,@(if (empty? atoms) '() `((atoms ,@atoms)))
                 ,@(if (empty? bonds) '() `((bonds ,@bonds))))))

  (define (sxml-from-atom a)
    `(atom (@ (label ,(interface:node.label a))
              (type ,(interface:atom.type a))
              (name ,(interface:atom.name a))
              ,@(if (= 1 (interface:atom.nsites a))
                    '()
                    `((nsites ,(number->string (interface:atom.nsites a))))))))

  (define (sxml-from-bond b)
    `(bond (@ (atoms ,(string-join (list (vector-ref b 0) (vector-ref b 1))
                                   " "))
              (order ,(vector-ref b 2)))))

  (define (sxml-from-transformation tr)
    (match-let*-values ([((cons t r)) tr]
                        [(t-sh t-data) (array->text t)]
                        [(r-sh r-data) (array->text r)])
      `(transformation (rotation ,r-data)
                       (translation ,t-data))))
  
  `(universe (@ ,@(if xml-id `((id ,xml-id)) '())
                (convention ,(interface:universe.convention universe))
                (cell_shape ,(interface:universe.cell-shape universe)))
             (symmetry_transformations
                 ,@(map sxml-from-transformation
                        (interface:universe.symmetry-transformations universe)))
             (molecules ,@(map sxml-from-molecule
                               (interface:universe.molecules universe)))))

(define (sxml-from-configuration configuration xml-id refs)
  (let*-values ([(cp-sh cp-data) (array->text (interface:configuration.cell-parameters configuration))]
                [(pos-sh pos-data) (array->text (interface:configuration.positions configuration))]
                [(ntype) (interface:configuration.number-type configuration)]
                [(universe) (interface:configuration.universe configuration)]
                [(universe-sxml) (if (dict-has-key? refs universe)
                                     `(universe (@ (ref ,(dict-ref refs universe))))
                                     (sxml-from-universe universe #f refs))])
    `(configuration (@ ,@(if xml-id `((id ,xml-id)) '()))
                    ,universe-sxml
                    (cell_parameters (@ (shape ,cp-sh)) ,cp-data)
                    (positions (@ (type ,ntype)) ,pos-data))))


(define (mosaic->sxml xml-id item [refs #hash()])
  (cond
   [(interface:universe? item) (sxml-from-universe item xml-id refs)]
   [(interface:configuration? item) (sxml-from-configuration item xml-id refs)]
   [else (error (format "unknown item type: ~a" item))]))


(define (mosaic-sequence->sxml-sequence seq)
  (in-generator
   (let ([refs (make-hash)])
     (~for ([($: xml-id item) seq])
       (let ([sxml (mosaic->sxml xml-id item refs)])
         (hash-set! refs item xml-id)
         (yield sxml))))))

(define (items-to-xml item-seq port)
  (displayln "<?xml version=\"1.0\" encoding=\"utf-8\"?>" port)
  (display (format "<mosaic version=\"~a.~a\">"
                   interface:mosaic-version-major
                   interface:mosaic-version-minor)
           port)
  (for ([sxml-item (mosaic-sequence->sxml-sequence item-seq)])
    (srl:sxml->xml-noindent sxml-item port))
    (displayln "</mosaic>" port))
