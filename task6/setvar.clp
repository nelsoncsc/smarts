(defrule set-variables
(load)
=>
(bind ?PO (read data))
(assert (PO ?PO))
(bind ?OR (read data))
(assert (OR ?OR))
(bind ?damping (read data))
(assert (damping ?damping))
(bind ?T (read data))
(assert (T ?T))
(bind ?Tr (read data))
(assert (Tr ?Tr)))