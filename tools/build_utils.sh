#!/bin/bashpp

ME=$(basename "$0")

if [[ $CI ]]; then
	TR=travis_retry
	#__QUIET = --quiet
fi

function download_boost {
	echo "$ME:$FUNCNAME: BEGIN in $(pwd)"

	if [[ $1 ]]; then __DIR=$1; else echo "$ME:$FUNCNAME: \$1 must be location"; return; fi;
	if [[ $2 ]]; then __VERSION=$2; else __VERSION=1.58.0; fi

	if [[ $__VERSION ]]; then
		BOOST_DIR=${__DIR}/boost-${__VERSION}

		if [[ ! -f $BOOST_DIR/.is_downloaded ]]; then
			mkdir -p "$BOOST_DIR"
			if [[ $__VERSION == "trunk" ]]; then
				BOOST_URL="http://github.com/boostorg/boost.git"
				$TR git clone --depth 1 --recursive "$__QUIET" $BOOST_URL "$BOOST_DIR" || exit 1
				(cd "$BOOST_DIR" && ./bootstrap.sh && ./b2 headers)
			else
				BOOST_URL="http://sourceforge.net/projects/boost/files/boost/${__VERSION}/boost_${__VERSION//\./_}.tar.gz"
				mkdir -p "${BOOST_DIR}"
				echo "$ME:$FUNCNAME: DOWNLOADING $BOOST_URL"
          		{ travis_retry wget ${__QUIET} -O - ${BOOST_URL} | tar --strip-components=1 -xz -C ${BOOST_DIR}; } || exit 1
          		echo "$ME:$FUNCNAME: DONE DOWNLOADING $BOOST_URL"
			fi
			touch "$BOOST_DIR/.is_downloaded"
		else
			echo "$ME:$FUNCNAME: ${BOOST_DIR}/.is_downloaded exists, skipping download"
		fi
	else
		echo "$ME:$FUNCNAME: no \$__VERSION (arg 2), doing nothing"
	fi
	echo "$ME:$FUNCNAME: END in $(pwd)"
}

function get_boost {
	echo "$ME:$FUNCNAME: BEGIN in $(pwd)"
	if [[ $1 ]]; then __DIR=$1; else echo "$ME:$FUNCNAME: $1 must be location"; return; fi;
	if [[ $2 ]]; then __VERSION=$2; else __VERSION=1.58.0; fi

	if [[ $__VERSION ]]; then
		if [[ ! -f $BOOST_DIR/.is_built ]]; then
			download_boost "$__DIR" $__VERSION
			(cd "$BOOST_DIR" && ./bootstrap.sh && ./b2 -d0 --prefix=./ --shared \
				--with-system --with-iostreams install && touch .is_built)
		else
			echo "$ME:$FUNCNAME: ${BOOST_DIR}/.is_built exists, skipping build"
		fi
		export CMAKE_OPTIONS=" -DBOOST_ROOT=${BOOST_DIR}"
	else
		echo "$ME:$FUNCNAME: no \$__VERSION (arg 2), doing nothing"
	fi
	echo "$ME:$FUNCNAME: END in $(pwd)"
}

function get_cmake {
	echo "$ME:$FUNCNAME: BEGIN in $(pwd)"
	if [[ $1 ]]; then __DIR=$1; else echo "$ME:$FUNCNAME: \$1 must be location"; return; fi;
	if [[ ! -f $__DIR/cmake/.is_downloaded ]]; then
		CMAKE_URL="https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.tar.gz"
		echo "$ME:$FUNCNAME: DOWNLOADING $CMAKE_URL"
		mkdir -p "$__DIR/cmake" && $TR wget --no-check-certificate "$__QUIET" -O - ${CMAKE_URL} \
			| tar --strip-components=1 -xz -C "$__DIR/cmake"
		echo "$ME:$FUNCNAME: DONE DOWNLOADING $CMAKE_URL"
		export PATH=${__DIR}/cmake/bin:${PATH}
		touch "$__DIR/cmake/.is_downloaded"
	else
		echo "$ME:$FUNCNAME: ${__DIR}/cmake/.is_downloaded exists, skipping download"
	fi
	echo "$ME:$FUNCNAME: END in $(pwd)"
}

function travis_get_cmake {
	echo "$ME:$FUNCNAME: BEGIN in $(pwd)"
	if [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
		get_cmake "$1" "$2"
	elif [ "$TRAVIS_OS_NAME" = 'osx' ]; then
		brew install cmake
		brew upgrade cmake
		export PATH=${PATH}:/Users/travis/Library/Python/2.7/bin
	else
		echo "$ME:$FUNCNAME: not on osx or linux???"
	fi
	echo "$ME:$FUNCNAME: END in $(pwd)"
}

function get_clang {
	echo "$ME:$FUNCNAME: BEGIN in $(pwd)"
	if [[ $1 ]]; then	__DIR=$1; else echo "$ME:$FUNCNAME: \$1 must be location"; return; fi;
	if [[ $2 ]]; then	__VERSION=$2; else __VERSION=3.9.0; fi

	if [[ $__VERSION ]]; then
		LLVM_DIR=${__DIR}/llvm-${__VERSION}
		if [[ ! -f $LLVM_DIR/.is_built ]]; then
			mkdir -p "$LLVM_DIR"
			LLVM_URL="http://llvm.org/releases/${__VERSION}/llvm-${__VERSION}.src.tar.xz"
			LIBCXX_URL="http://llvm.org/releases/${__VERSION}/libcxx-${__VERSION}.src.tar.xz"
			LIBCXXABI_URL="http://llvm.org/releases/${__VERSION}/libcxxabi-${__VERSION}.src.tar.xz"
			CLANG_URL="http://llvm.org/releases/${__VERSION}/clang+llvm-${__VERSION}-x86_64-linux-gnu-ubuntu-14.04.tar.xz"

			mkdir -p "${LLVM_DIR}" "${LLVM_DIR}/build" "${LLVM_DIR}/projects/libcxx" \
				"${LLVM_DIR}/projects/libcxxabi" "${LLVM_DIR}/clang"
			if [[ ! -f $LLVM_DIR/.is_downloaded ]]; then
			  $TR wget --quiet -O - ${LLVM_URL}      | tar --strip-components=1 -xJ -C ${LLVM_DIR}
   	    $TR wget --quiet -O - ${LIBCXX_URL}    | tar --strip-components=1 -xJ -C ${LLVM_DIR}/projects/libcxx
       	$TR wget --quiet -O - ${LIBCXXABI_URL} | tar --strip-components=1 -xJ -C ${LLVM_DIR}/projects/libcxxabi
       	$TR wget --quiet -O - ${CLANG_URL}     | tar --strip-components=1 -xJ -C ${LLVM_DIR}/clang
				touch "$LLVM_DIR/.is_downloaded"
			else
				echo "$ME:$FUNCNAME: ${LLVM_DIR}/.is_downloaded exists, skipping llvm download"
      fi
			(cd "${LLVM_DIR}/build" && cmake .. -DCMAKE_INSTALL_PREFIX="${LLVM_DIR}/install" -DCMAKE_CXX_COMPILER=clang++)
			(cd "${LLVM_DIR}/build/projects/libcxx" && make install -j2)
			(cd "${LLVM_DIR}/build/projects/libcxxabi" && make install -j2)
			touch "$LLVM_DIR/.is_built"
		else
			echo "$ME:$FUNCNAME: ${LLVM_DIR}/.is_built exists, skipping llvm build"
		fi
	else
		echo "$ME:$FUNCNAME: no \$__VERSION (arg2), doing nothing"
	fi
	export CC=clang
	export CXX=clang++
	export CXXFLAGS="-nostdinc++ -isystem ${LLVM_DIR}/install/include/c++/v1"
	export LDFLAGS="-L${LLVM_DIR}/install/lib -lc++ -lc++abi"
	export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${LLVM_DIR}/install/lib"
	export PATH="${LLVM_DIR}/clang/bin:${PATH}"
	echo "$ME:$FUNCNAME: END in $(pwd)"
}


function get_doxygen {
	echo "$ME:$FUNCNAME: BEGIN in $(pwd)"
	DOXYGEN_URL="http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.11.linux.bin.tar.gz"
	echo "$ME:$FUNCNAME: DOWNLOADING $DOXYGEN_URL"
	mkdir doxygen && $TR wget "$__QUIET" -O - ${DOXYGEN_URL} | tar --strip-components=1 -xz -C doxygen
	echo "$ME:$FUNCNAME: DONE DOWNLOADING $DOXYGEN_URL"
	echo "$ME:$FUNCNAME: add doxygen to path"
	export PATH=${__DIR}/doxygen/bin:${PATH}
	doxygen --version
	echo "$ME:$FUNCNAME: END in $(pwd)"
}
