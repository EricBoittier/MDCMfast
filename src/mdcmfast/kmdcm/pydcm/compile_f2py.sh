python -m numpy.f2py -c -m dcm_fortran dcm_fortran.F90 --fcompiler=gnu95 \
  -I$HOME/.local/share/uv/python/cpython-3.12.3-linux-x86_64-gnu/include
 #--debug-capi

