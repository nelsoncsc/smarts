(defrule get-kc-small
(Kc-is-small)
=>
(printout Kc 0 crlf))

(defrule get-kc-big
(Kc-is-big)
=>
(printout Kc 1 crlf))

(defrule get-ti-small
(Ti-is-small)
=>
(printout Ti 0 crlf))

(defrule get-ti-big
(Ti-is-big)
=>
(printout Ti 1 crlf))

(defrule get-td-small
(Td-is-small)
=>
(printout Td 0 crlf))

(defrule get-td-big
(Td-is-big)
=>
(printout Td 1 crlf))


