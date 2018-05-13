all : index.html

index.html : LUAD_project.knit.md
	Rscript -e "rmarkdown::render('$<')"

LUAD_project.knit.md : LUAD_project.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

