ADD_EXECUTABLE(jump_diffusion jump_diffusion.cpp)
ADD_EXECUTABLE(tutorial tutorial.cpp)
ADD_EXECUTABLE(unequal_borrowing_lending_rates unequal_borrowing_lending_rates.cpp)
ADD_EXECUTABLE(vanilla_options vanilla_options.cpp)

TARGET_LINK_LIBRARIES(jump_diffusion ${JSONCPP_LIBRARY})
TARGET_LINK_LIBRARIES(unequal_borrowing_lending_rates ${JSONCPP_LIBRARY})
TARGET_LINK_LIBRARIES(vanilla_options ${JSONCPP_LIBRARY})
