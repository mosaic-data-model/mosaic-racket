#lang racket

(require rackunit
         math/array
         generic-bind
         "../xml.rkt"
         "../validation.rkt"
         (prefix-in interface: "../interface.rkt")
         (prefix-in model: "../model.rkt"))

(provide xml-tests)

(define xml-test-files '("ala_dipeptide.xml" "water.xml" "3ZQ8.xml"))

(define (check-items items)
  (check-equal?
   (filter (λ (m) (not (equal? (cdr m) "unknown item type")))
           (messages (validate-items (sequence-map cdr items)) '()))
   '()))

(define-simple-check (check-equivalent? item1 item2)
  (empty? (interface:diffs item1 item2)))

(define (wrap-sxml-list sxml-list)
  (let ([refs (make-hash)])
    (for/list ([node sxml-list])
      (match-let ([(cons xml-id item) (sxml->mosaic node refs)])
        (hash-set! refs xml-id item)
        item))))

(define xml-tests
  (test-suite "xml"

(test-case "items"
  (for ([filename xml-test-files])
    (check-items (items-from-xml (open-input-file filename)))))

(test-case "conversion"
  (for ([filename xml-test-files])
    (let* (; Read XML file and return SXML nodes with Mosaic tags
           [sxml-items-with-ids
            (filter (λ (x)
                      (or (interface:universe? (cdr x))
                          (interface:configuration? (cdr x))))
                    (call-with-input-file filename
                      (λ (in) (sequence->list (items-from-xml in)))))]
           ; Create Mosaic model data items from SXML
           [model-items-with-ids (map (λ (x)
                                        (cons (car x)
                                              (model:make-data-item (cdr x))))
                                      sxml-items-with-ids)]
           ; Create Mosaic SXML items from Mosaic model items
           [test-items (wrap-sxml-list
                        (sequence->list
                         (mosaic-sequence->sxml-sequence model-items-with-ids)))]
           ; Remove xml-ids
           [sxml-items (map cdr sxml-items-with-ids)]
           [model-items (map cdr model-items-with-ids)]
           ; Create XML string from item sequence
           [xml-string (let ([out (open-output-string)])
                         (items-to-xml model-items-with-ids out)
                         (get-output-string out))]
           ; Create item sequence from XML string
           [recovered-items (map cdr (sequence->list
                                      (items-from-xml
                                       (open-input-string xml-string))))])
      (for ([sxml sxml-items]
            [model model-items]
            [test test-items]
            [recovered recovered-items])
        (check-equivalent? sxml model)
        (check-equivalent? sxml test)
        (check-equivalent? sxml recovered)
        (check-equivalent? model test)
        (check-equivalent? model recovered)
        (check-equivalent? test recovered)))))))
