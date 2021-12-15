GIT_SHA1:=$(shell git rev-parse HEAD)
MACRO_DEFS+=-DGIT_SHA1='"$(GIT_SHA1)"'

# Version string, including current commit if not on a release, plus
# '-dirty' if there are uncommitted changes
GIT_VERSION := $(shell git describe --tags --long --dirty --match "[0-9]*.[0-9]*.[0-9]*")
MACRO_DEFS += -DGIT_VERSION='"$(GIT_VERSION)"'

# Find if there are any modified tracked files (except Makefile.depend)
ifeq ($(shell git status --short -uno -- . | wc -l), 0)
	GIT_STATE:="clean"
	MACRO_DEFS+=-DGIT_STATE='$(GIT_STATE)'
$(shell touch $(PATCHFILE))
else
	GIT_STATE:="modified"
	MACRO_DEFS+=-DGIT_STATE='$(GIT_STATE)'
