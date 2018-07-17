(TeX-add-style-hook
 "introduction"
 (lambda ()
   (TeX-run-style-hooks
    "introduction/motivation"
    "introduction/background"
    "introduction/contributions"))
 :latex)

