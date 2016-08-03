(defrule openFile
(load)
=>
(open "data.csv" data "r")
(open "Kc.txt" Kc "w")
(open "Ti.txt" Ti "w")
(open "Td.txt" Td "w"))

