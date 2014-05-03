all:
	$(MAKE) -C tests

doc:
	doxygen

clean:
	$(RM) -r doc
	$(MAKE) clean -C tests

.PHONY: all doc clean

