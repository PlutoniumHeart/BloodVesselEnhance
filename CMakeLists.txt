cmake_minimum_required( VERSION 2.8 )

project ( BloodVesselEnhance )

# SET(CMAKE_CXX_FLAGS "-std=c++11")

find_package(CUDA REQUIRED)
if(CUDA_FOUND)
  cuda_include_directories(${CUDA_SDK_ROOT_DIR}/common/inc ${CUDA_TOOLKIT_INCLUDE} ${CUDA_INCLUDE_DIRS})
  include_directories(${CUDA_SDK_ROOT_DIR}/common/inc)
  link_directories(${CUDA_SDK_ROOT_DIR}/common/lib/x64)
  link_directories( ${CUDA_SDK_ROOT_DIR})
  link_directories( ${CUDA_TOOLKIT_ROOT_DIR}/lib/x64 )
  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
else(CUDA_FOUND)
  message(FATAL_ERROR "CUDA not found. Please set the related directories manually.")
endif(CUDA_FOUND)

find_package(ITK REQUIRED)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
else(ITK_FOUND)
  message(FATAL_ERROR "ITK not found. Please set ITK_DIR. ")
endif(ITK_FOUND)

include_directories(src/include)

set(SOURCES src/main.cpp src/include/BaseImage.h src/BaseImage.cpp
            src/include/BaseImageSeries.h src/BaseImageSeries.cpp
            src/include/DicomSlice.h src/DicomSlice.cpp
            src/include/DicomSeries.h src/DicomSeries.cpp
            src/include/PngImage.h src/PngImage.cpp
            src/include/PngSeries.h src/PngSeries.cpp
            src/include/JpegImage.h src/JpegImage.cpp
            src/include/JpegSeries.h src/JpegSeries.cpp
            src/include/StructureClassify.h src/StructureClassify.cpp
            src/include/Matrix.h src/Matrix.cpp
            src/include/BaseFilter.h src/BaseFilter.cpp
            src/include/RecursiveGaussian.h src/RecursiveGaussian.cpp
            src/include/GPUMemOps.cuh src/GPUMemOps.cu
            src/cudaRecursiveGaussian_kernel.cu
)

set(ALL_LIBS ${ITK_LIBRARIES} ${CUDA_LIBRARIES})

cuda_add_executable(BloodVesselEnhance ${SOURCES})

target_link_libraries(BloodVesselEnhance ${ALL_LIBS})
