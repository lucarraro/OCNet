# Makefile to be placed on the level of OCNet
# simplifies structure


# the following could be put in file Makefile.def
# and `include Makefile.def` here
SHELL := /bin/bash
UNAME = $(shell uname -n)
ifeq ($(UNAME),illnau)
   R := /usr/lib/R-devel/bin/R
   RSCRIPT := /usr/lib/R-devel/bin/Rscript
else
  R :=  R
  RSCRIPT :=  Rscript
endif



VERSION := 0.4
 # "$(shell grep 'Version' OCNet/DESCRIPTION)"
TARPACK = OCNet_$(VERSION).tar.gz

CHECKARGS =   --as-cran # --use-gct --use-valgrind
# --run-dontrun     do run \dontrun sections in the Rd files
# --run-donttest    do run \donttest sections in the Rd files
# --use-gctuse      'gctorture(TRUE)' when running examples/tests

##  manual/automatic change of version numbers and date/timestamps
update:
	if [ "Version: $(VERSION)" != "$(shell grep 'Version:' OCNet/DESCRIPTION)" ] ; then \
	echo "Missmatched Version number!"; \
	fi

## The following could be used to automatically set the date in Description or in the vignette. 
#	if [ "Date: $(shell date +'%F')" != "$(shell grep 'Date' OCNet/DESCRIPTION)" ] ; then \
#	cd OCNet && sed -i -r -- 's/^Date:.*/Date: '`date +'%F'`'/g' DESCRIPTION ; \
#	fi


tar: tar/$(TARPACK)
tar/$(TARPACK): $(shell find OCNet -type f) Makefile
	$(MAKE) update
	mkdir -p tar
	cd tar && $(R) CMD build ../OCNet

## install to 'lib' ------------------------------------
lib: 
	mkdir -p lib
	$(R) CMD INSTALL -l lib OCNet

lib-systemOCNet:
	$(RSCRIPT) -e "install.packages(\"tar/OCNet_$(VERSION).tar.gz\", repos=NULL)"


check: tar lib-systemOCNet
	cd tar && $(R) CMD check  $(CHECKARGS) $(TARPACK)

## winbuilder
check-win: update
	$(RSCRIPT) -e "devtools::check_win_release(pkg = \"OCNet\")"

check-win-old: update
	$(RSCRIPT) -e "devtools::check_win_oldrelease(pkg = \"OCNet\")"

check-win-devel: update
	$(RSCRIPT) -e "devtools::check_win_devel(pkg = \"OCNet\")"


## test package -------------------------------------------
## run tests including skip_on_cran() tests
test-package: lib
	$(RSCRIPT) -e "library(\"OCNet\", lib.loc = \"lib\"); devtools::test(\"OCNet\")"

## run examples including "dontrun example"
test-examples: lib
	$(RSCRIPT) -e "library(\"OCNet\", lib.loc = \"lib\"); devtools::run_examples(\"OCNet\", run = FALSE)"
	rm -f Rplots*.pdf



finalizer: 
	cd tar && $(R) CMD check $(TARPACK)
	cd tar && cp OCNet.Rcheck/tests/test-package.Rout  ../OCNet/tests/test-package.Rout.save
	cd tar && cp OCNet.Rcheck/OCNet-Ex.Rout          ../OCNet/tests/Examples/OCNet-Ex.Rout.save
	cd tar && $(R) CMD check $(TARPACK)

clean:
	rm -f  Rplots*.pdf .RData .Rhistory .RData


.PHONY:	update \
	tar lib lib-systemOCNet \
	check check-win check-win-old  check-win-devel\
	test-package test-examples \
        clean \
	finalizer
