include(UseLATEX)


ADD_LATEX_DOCUMENT(userManual.tex
   INPUTS figs/TemperatureProfiles.tex
   BIBFILES userManual.bib
   CONFIGURE userManual.tex
   IMAGE_DIRS figs
   DEFAULT_PDF
)

