#!/bin/bash -e

show_help() {
    echo "Usage: source $(basename "$BASH_SOURCE") [OPTION]"
    echo "Valid options:"
    echo "  --heffte_build_type=<fftw|cufft|rocfft>"
    echo "  --kokkos_build_type=<serial|cuda|hip>"
    echo "  --help: Display this help message"
    return 1
}

# Initialize variables with default values
heffte_build_type=""
kokkos_build_type=""

# Define arrays of valid options
valid_heffte_build_types=("fftw" "cufft" "rocfft")
valid_kokkos_build_types=("serial" "cuda" "hip")

# Parse command line arguments
for arg in "$@"; do
    case "$arg" in
        --heffte_build_type=*)
            option="${arg#*=}"
            if [[ " ${valid_heffte_build_types[*]} " == *" $option "* ]]; then
                heffte_build_type="$option"
            else
                echo "Error: Invalid --heffte_build_type specified."
                show_help
                return 1
            fi
            ;;
        --kokkos_build_type=*)
            option="${arg#*=}"
            if [[ " ${valid_kokkos_build_types[*]} " == *" $option "* ]]; then
                kokkos_build_type="$option"
            else
                echo "Error: Invalid --kokkos_build_type specified."
                show_help
                return 1
            fi
            ;;
        --help)
            show_help
            return 1
            ;;
        *)
            echo "Error: Invalid argument or value specified."
            show_help
            return 1
            ;;
    esac
done

# Check if required options are specified
if [ -z "$heffte_build_type" ] || [ -z "$kokkos_build_type" ]; then
    echo "Error: --heffte_build_type and --kokkos_build_type are required options."
    show_help
    return 1
fi

# If both arguments are valid, you can use them in your script as needed
echo "Heffte Build Type: $heffte_build_type"
echo "Kokkos Build Type: $kokkos_build_type"

# Determine the directory of the current script
SCRIPT_DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
echo "Script directory: ${SCRIPT_DIR}"

# Determine the parent directory of the script's directory
PARENT_DIR=$(dirname $(dirname "${SCRIPT_DIR}"))

# --------building heffte
HEFFTE_CONFIG_SCRIPT="$PARENT_DIR/scripts/install-scripts/install_heffte.sh"
source "$HEFFTE_CONFIG_SCRIPT" --$heffte_build_type

# --------building kokkos
KOKKOS_CONFIG_SCRIPT="$PARENT_DIR/scripts/install-scripts/install_kokkos.sh"
source "$KOKKOS_CONFIG_SCRIPT" --$kokkos_build_type

# --------building hdf5
HDF5_CONFIG_SCRIPT="$PARENT_DIR/scripts/install-scripts/install_hdf5.sh"
source "$HDF5_CONFIG_SCRIPT"

# Check if the 'MATAR' directory exists in the parent directory; if not, clone it
MATAR_DIR="$PARENT_DIR/MATAR"
if [ ! -d "$MATAR_DIR" ]; then
  echo "Directory 'MATAR' does not exist in '$PARENT_DIR', downloading 'MATAR'...."
  git clone https://github.com/lanl/MATAR.git "$MATAR_DIR"
else
  echo "Directory 'MATAR' exists in '$PARENT_DIR', skipping 'MATAR' download"
fi

# --------building EVPFFT
EVPFFT_SOURCE_DIR="$PARENT_DIR/src"
EVPFFT_BUILD_DIR="$PARENT_DIR/evpfft_heffte_${heffte_build_type}_kokkos_${kokkos_build_type}"
HEFFTE_INSTALL_DIR="$PARENT_DIR/heffte/install_heffte_$heffte_build_type"
KOKKOS_INSTALL_DIR="$PARENT_DIR/kokkos/install_kokkos_$kokkos_build_type"
HDF5_INSTALL_DIR="$PARENT_DIR/hdf5/install"
MATAR_SOURCE_DIR="$PARENT_DIR/MATAR/src"

# Configure EVPFFT using CMake
cmake_options=(
    -D CMAKE_BUILD_TYPE=Release
    -D CMAKE_PREFIX_PATH="$HEFFTE_INSTALL_DIR;$KOKKOS_INSTALL_DIR;$HDF5_INSTALL_DIR"
    -D CMAKE_CXX_FLAGS="-I$MATAR_SOURCE_DIR"
    -D ENABLE_PROFILING=ON
)

if [ "$heffte_build_type" == "fftw" ]; then
    cmake_options+=(
        -D USE_FFTW=ON
    )   
elif [ "$heffte_build_type" == "cufft" ]; then
    cmake_options+=(
        -D USE_CUFFT=ON
    )   
elif [ "$heffte_build_type" == "rocfft" ]; then
    cmake_options+=(
      -D USE_ROCFFT=ON
    )   
fi

# Configure EVPFFT
cmake "${cmake_options[@]}" -B "$EVPFFT_BUILD_DIR" -S "$EVPFFT_SOURCE_DIR"

# Build kokkos
echo "Building EVPFFT..."
make -C "$EVPFFT_BUILD_DIR" -j
