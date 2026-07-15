install:
	Rscript -e "devtools::install()"

check:
	Rscript -e "devtools::check()"

build:
	Rscript -e "pkgdown::build_site()"

