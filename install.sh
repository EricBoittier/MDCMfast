cd src/mdcmfast/bin/
make 

cd ../kmdcm/pydcm


# 1. Make sure the expected include path exists
mkdir -p ~/.local/include
ln -s ~/.local/share/uv/python/cpython-3.12.3-linux-x86_64-gnu/include/python3.12 ~/.local/include/python3.12 2>/dev/null || true

# 2. Export compiler flags to find Python.h
export CFLAGS="-I$HOME/.local/include/python3.12"
export CPPFLAGS="$CFLAGS"

# 3. (Optional) Clean and reattempt the build
rm -rf bbdir  # Clear previous failed build

bash compile_f2py.sh


