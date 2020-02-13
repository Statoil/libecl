#!/bin/bash

set -e
set -x

PROJECT=libecl
RELEASE_PATH=${KOMODO_ROOT}/${RELEASE_NAME}
echo "unit testing ${PROJECT} against ${RELEASE_NAME}"

source ${RELEASE_PATH}/enable
source ${DEVTOOL}/enable
GCC_VERSION=7.3.0 CMAKE_VERSION=3.10.2 source ${SDPSOFT}/env.sh

GIT=${SDPSOFT}/bin/git

EV=$(cat ${RELEASE_PATH}/${RELEASE_NAME} | grep "${PROJECT}:" -A2 | grep "version:")
EV=($EV)    # split the string "version: vX.X.X"
EV=${EV[1]} # extract the version
echo "Using ${PROJECT} version ${EV}"
$GIT checkout $EV

rm -rf build
mkdir build
pushd build
cmake .. -DBUILD_TESTS=ON \
-DENABLE_PYTHON=ON \
-DBUILD_APPLICATIONS=ON \
-DINSTALL_ERT_LEGACY=ON \
-DERT_USE_OPENMP=ON \
-DCMAKE_INSTALL_PREFIX=install \
-DCMAKE_C_FLAGS='-Werror=all' \
-DCMAKE_CXX_FLAGS='-Wno-unused-result'
make -j 12
#removing built libs in order to ensure we are using libs from komodo
rm -r lib64

echo "running ctest"
ctest --output-on-failure
popd
