
;; Overshoot Percentual Range Small
(defrule PO-small
(PO ?PO)
(test (< ?PO 50))
 =>
(assert (PO-is-small)))

;; Overshoot Percentual Range Big
(defrule PO-big
(PO ?PO)
(test (>= ?PO 50))
 =>
(assert (PO-is-big)))

;; Overshoot Ratio Range Small
(defrule OR-small
(OR ?OR)
(test (< ?OR 0.5))
 =>
(assert (OR-is-small)))

;; Overshoot Ratio Range Big
(defrule OR-big
(OR ?OR)
(test (>= ?OR 0.5))
 =>
(assert (OR-is-big)))

;; Damping Range Small
(defrule damping-small
(damping ?damping)
(test (< ?damping 0.5))
 =>
(assert (damping-is-small)))

;; Damping Range Big
(defrule damping-big
(damping ?damping)
(test (>= ?damping 0.5))
 =>
(assert (damping-is-big)))

;; Period Range Small
(defrule T-small
(T ?T)
(test (< ?T 50))
 =>
(assert (T-is-small)))

;; Period Range Big
(defrule T-big
(T ?T)
(test (>= ?T 50))
 =>
(assert (T-is-big)))

;; Rise Time Range Small
(defrule Tr-small
(Tr ?Tr)
(test (< ?Tr 50))
 =>
(assert (Tr-is-small)))

;; Rise Time Range Big
(defrule Tr-big
(Tr ?Tr)
(test (>= ?Tr 50))
 =>
(assert (Tr-is-big)))

;; If PO is Small then Kc is Small
(defrule PO-Kc-Small

(PO-is-small)
=> 
(assert (Kc-is-small)))

;; If PO is Big then Kc is Big
(defrule PO-Kc-Big

(PO-is-big)
=> 
(assert (Kc-is-big)))

;; If damping is Small then Kc is Small
(defrule damping-Kc-Small

(damping-is-small)
=> 
(assert (Kc-is-small)))

;; If damping is Big then Kc is Big
(defrule damping-Kc-Big

(damping-is-big)
=> 
(assert (Kc-is-big)))

;; If T is Big then Kc is Small
(defrule T-Kc-Small

(T-is-big)
=> 
(assert (Kc-is-small)))

;; If T is Small then Kc is Big
(defrule T-Kc-Big

(T-is-small)
=> 
(assert (Kc-is-big)))

;; If PO is Big then Ti is Small
(defrule PO-Ti-Small

(PO-is-big)
=> 
(assert (Ti-is-small)))

;; If PO is Small then Ti is Big
(defrule PO-Ti-Big

(PO-is-small)
=> 
(assert (Ti-is-big)))

;; If damping is Big then Ti is Small
(defrule damping-Ti-Small

(damping-is-big)
=> 
(assert (Ti-is-small)))

;; If damping is Small then Ti is Big
(defrule damping-Ti-Big

(damping-is-small)
=> 
(assert (Ti-is-big)))

;; If OR is Big then Ti is Small
(defrule OR-Ti-Small

(OR-is-big)
=> 
(assert (Ti-is-small)))

;; If OR is Small then Ti is Big
(defrule OR-Ti-Big

(OR-is-small)
=> 
(assert (Ti-is-big)))

;; If Tr is Big then Td is Small
(defrule Tr-Td-Small

(Tr-is-big)
=> 
(assert (Td-is-small)))

;; If Tr is Small then Td is Big
(defrule Tr-Td-Big

(Tr-is-small)
=> 
(assert (Td-is-big)))
