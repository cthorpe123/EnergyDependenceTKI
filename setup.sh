source /grid/fermiapp/products/uboone/setup_uboone.sh
setup uboonecode v08_00_00_69 -q e17:prof
export ROOT_INCLUDE_PATH=$PWD/Include:$ROOT_INCLUDE_PATH
export LD_LIBRARY_PATH=$PWD/Include:$LD_LIBRARY_PATH
