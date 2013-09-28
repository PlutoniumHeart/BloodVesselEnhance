cmake_minimum_required( VERSION 2.8 )

project ( Opti511_Final )

find_package(ITK REQUIRED)
if(ITK_FOUND)
    include(${ITK_USE_FILE})
else(ITK_FOUND)
    message(FATAL_ERROR "ITK not found. Please set ITK_DIR. ")
endif(ITK_FOUND)

find_package(CUDA REQUIRED)
if(CUDA_FOUND)
	cuda_include_directories(${CUDA_SDK_ROOT_DIR}/common/inc ${CUDA_TOOLKIT_INCLUDE} ${CUDA_INCLUDE_DIRS})
	include_directories(${CUDA_SDK_ROOT_DIR}/common/inc)
	link_directories(${CUDA_SDK_ROOT_DIR}/common/lib/x64)
	set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
else(CUDA_FOUND)
	message(FATAL_ERROR "CUDA not found. Please set the related directories manually.")
endif(CUDA_FOUND)

set(SOURCES src/main.cpp src/BaseImage.h src/BaseImage.cpp
            src/DicomSlice.h src/DicomSlice.cpp
            src/PngImage.h src/PngImage.cpp
            src/JpegImage.h src/JpegImage.cpp
)

set(ALL_LIBS ${ITK_LIBRARIES} ${CUDA_LIBRARIES})

cuda_add_executable(Opti511_Final ${SOURCES})

target_link_libraries(Opti511_Final ${ALL_LIBS})