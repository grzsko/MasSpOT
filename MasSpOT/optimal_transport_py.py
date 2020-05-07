import os
from cffi import FFI
ffi = FFI()
ffi.cdef("""double calc_distance_c(double*, double*, int, double*, double*, int,
          double, double, double, double, double);""")

base_dir = os.path.dirname(os.path.abspath(__file__)),
base_dir = os.path.join(base_dir[0], "..")

for file in os.listdir(base_dir):
    print(file)
    if file.startswith("MasSpOTCppToPy") and file.endswith(".so"):
        so_path = os.path.join(base_dir, file)
lib = ffi.dlopen(so_path)

# Above searches for .so files in installed package.


def spectral_distance(spec1, spec2, lam=1, eps=0.1, tol=1e-05, threshold=1e+02,
                   max_iter=500):

    mzs1, ints1 = zip(*[x for x in spec1.confs if x[1] > 0])
    mzs2, ints2 = zip(*[x for x in spec2.confs if x[1] > 0])

    ffi_mzs1 = ffi.new("double[]", mzs1)
    ffi_mzs2 = ffi.new("double[]", mzs2)
    ffi_ints1 = ffi.new("double[]", ints1)
    ffi_ints2 = ffi.new("double[]", ints2)


    res = lib.calc_distance_c(
        ffi_mzs1, ffi_ints1, len(mzs1), ffi_mzs2, ffi_ints2, len(mzs2),
        lam, eps, tol, threshold, max_iter)

    return res

