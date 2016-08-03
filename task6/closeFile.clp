(defrule closeData
(close data)
=>
(close data))
(defrule closeKc
(= read Kc EOF)
=>
(close Kc))
(defrule closeTi
(close Ti)
=>
(close Ti))

(defrule closeTd
(close Td)
=>
(close Td))

