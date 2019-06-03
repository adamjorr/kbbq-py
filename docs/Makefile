# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
CLIHELP       = $(SOURCEDIR)/cli/help.rst

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile $(CLIHELP)
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

$(CLIHELP):
	kbbq -h > $(CLIHELP) &&\
	echo -e '\n\n' >> $(CLIHELP) &&\
	kbbq recalibrate -h >> $(CLIHELP) &&\
	echo -e '\n\n' >> $(CLIHELP) &&\
	kbbq benchmark -h >> $(CLIHELP) &&\
	echo -e '\n\n' >> $(CLIHELP) &&\
	kbbq plot -h >> $(CLIHELP) &&\
	echo -e '\n\n' >> $(CLIHELP)
