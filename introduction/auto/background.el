(TeX-add-style-hook
 "background"
 (lambda ()
   (TeX-run-style-hooks
    "introduction/background_computational_complexity"
    "introduction/background_complex_geometry"
    "introduction/background_petrov_galerkin")
   (LaTeX-add-labels
    "sec:intro_background"))
 :latex)

