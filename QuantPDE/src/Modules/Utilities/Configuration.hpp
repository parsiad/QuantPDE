#ifndef QUANT_PDE_MODULES_CONFIGURATION_HPP
#define QUANT_PDE_MODULES_CONFIGURATION_HPP

#include <json/json.h>

#include <array>    // std::array
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
		const auto tmp = (type) \
				configuration.get(key, defaultValue).asType();\
		configuration[key] = tmp; \
		return tmp; \
	}

QUANT_PDE_CONFIGURATION_GET(getInt, int, asInt)
QUANT_PDE_CONFIGURATION_GET(getReal, Real, asDouble)
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

class ConfigurationHelpers final {

public:

	ConfigurationHelpers() = delete;
	ConfigurationHelpers (const ConfigurationHelpers &) = delete;
	ConfigurationHelpers &operator=(const ConfigurationHelpers &) = delete;

	template <Index Dimension>
	static RectilinearGrid<Dimension> getGrid(
		Configuration &configuration,
		const std::string &key,
		const RectilinearGrid<Dimension> &defaultValue
	) {

		if(!configuration.isMember(key)) {
			for(Index i = 0; i < Dimension; ++i) {
				for(
					Index j = 0;
					j < defaultValue[i].size();
					++j
				) {
					configuration[key][i][j]
							= defaultValue[i][j];
				}
			}
			return defaultValue;
		}

		Axis axes[Dimension];
		for(Index i = 0; i < Dimension; ++i) {
			axes[i].length = configuration[key][i].size();
			axes[i].n = new Real[axes[i].length];
			for(Index j = 0; j < axes[i].length; ++j) {
				axes[i].n[j] = (Real) configuration[key][i][j]
						.asDouble();
			}
		}

		return RectilinearGrid<Dimension>(axes);

	}

};

template <Index Dimension>
inline RectilinearGrid<Dimension> getGrid(
	Configuration &configuration,
	const std::string &key,
	const RectilinearGrid<Dimension> &defaultValue
) {
	return ConfigurationHelpers::getGrid<Dimension>(
		configuration,
		key,
		defaultValue
	);
}

}

#endif

