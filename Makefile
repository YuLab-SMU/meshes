PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rd check clean

alldocs: rd readme

rd:
	Rscript -e 'roxygen2::roxygenise(".")'

readme:
	Rscript -e 'rmarkdown::render("README.Rmd")'

build:
	cd ..;\
	R CMD build $(PKGSRC)

build2:
	cd ..;\
	R CMD build --no-build-vignettes $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	Rscript -e 'rcmdcheck::rcmdcheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

check2: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz;\

bioccheck:
	cd ..;\
	Rscript -e 'BiocCheck::BiocCheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

rmoldrelease:
	git branch -D RELEASE_3_11

release:
	git checkout RELEASE_3_12;\
	git fetch --all


gitmaintain:
	git gc --auto;\
	git prune -v;\
	git fsck --full

update:
	git fetch --all;\
	git checkout master;\
	git merge upstream/master;\
	git merge origin/master

push: update
	git push upstream master;\
	git push origin master

biocinit:
	git remote add upstream git@git.bioconductor.org:packages/$(PKGNAME).git;\
	git fetch --all
