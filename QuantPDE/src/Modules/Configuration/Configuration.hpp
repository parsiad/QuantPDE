#ifndef QUANT_PDE_MODULES_CONFIGURATION_HPP
#define QUANT_PDE_MODULES_CONFIGURATION_HPP

#include <json/json.h>

#include <cstdlib>  // std::exit
#include <fstream>  // std::ifstream
#include <iostream> // std::cin
#include <string>   // std::string
#include <unistd.h> // std::getopt

namespace QuantPDE {

typedef Json::Value Configuration;

#define QUANT_PDE_CONFIGURATION_GET(name, type, asType) \
	type name( \
		Configuration &configuration, \
		const std::string &key, \
		type defaultValue \
	) { \
		const auto tmp = configuration.get(key, defaultValue).asType();\
		configuration[key] = tmp; \
		return tmp; \
	}

QUANT_PDE_CONFIGURATION_GET(getInt, int, asInt)
QUANT_PDE_CONFIGURATION_GET(getDouble, double, asDouble)
QUANT_PDE_CONFIGURATION_GET(getBool, bool, asBool)

#undef QUANT_PDE_CONFIGURATION_GET

Configuration getConfiguration(int argc, char **argv) {

	bool input = false;

	char c;
	while((c = getopt(argc, argv, "hi")) != -1) {
		switch(c) {
			case 'h':
				std::cerr
					<< argv[0]
					<< " [CONFIGURATION_FILE | -i]"
					<< std::endl
				;
				std::exit(2);
				break;
			case 'i':
				input = true;
				break;
		}
	}

	Configuration configuration;

	if(input) {
		std::cin >> configuration;
	} else if(argc >= 2) {
		std::ifstream(argv[1], std::ifstream::binary) >> configuration;
	}

	return configuration;

}

}

#endif

