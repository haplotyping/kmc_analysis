SUBDIRS := query analysis python

.PHONY: all clean $(SUBDIRS)

all clean:
	mkdir -p bin
	mkdir -p lib
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir $@; \
	done

$(SUBDIRS):
	$(MAKE) -C $@


  


