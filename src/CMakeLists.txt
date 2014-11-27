set(BAMREAD_SRCS bam-read.cpp pileup.cpp segmenter.cpp)
set(BAMREAD_HEADERS pileup.h segmenter.h)

set(BAMSPLIT_SRCS bam-split.cpp)
set(BAMSPLIT_HEADERS bam-split.h)


include_directories(${GAT_SOURCE_DIR}/bamtools/include)
link_directories(${GAT_SOURCE_DIR}/bamtools/lib)

set(LIBRARIES ${ZLIB_LIBRARIES} ${GAT_SOURCE_DIR}/bamtools/lib/libbamtools.a)

add_executable(bam-read ${BAMREAD_SRCS} ${BAMREAD_HEADERS})
add_executable(bam-split ${BAMSPLIT_SRCS} ${BAMSPLIT_HEADERS})

target_link_libraries(bam-read ${LIBRARIES})
target_link_libraries(bam-split ${LIBRARIES})
