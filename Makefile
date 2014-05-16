all:
	$(MAKE) -C Tests

Documentation:
	doxygen

clean:
	$(RM) -r Documentation
	$(MAKE) clean -C Tests

.PHONY: all clean

