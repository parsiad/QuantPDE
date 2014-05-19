all:
	$(MAKE) -C test

doc:
	doxygen

clean:
	$(RM) -r doc
	$(MAKE) clean -C test

.PHONY: all clean

