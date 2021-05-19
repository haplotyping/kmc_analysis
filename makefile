SUBDIRS := query analysis

.PHONY: all clean $(SUBDIRS)

all clean:
	mkdir -p bin
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir $@; \
	done

$(SUBDIRS):
	$(MAKE) -C $@


  


