
using Libdl

function main()
    R = 0.07
    E = 2.0
    B = 0.17
    D = 0.42
    G = 0.09
    H = 0.1
    Q = 0.4
    
    lib = dlopen("dynasty.dylib", RTLD_LAZY | RTLD_GLOBAL)
    func_ptr = dlsym(lib, "integrate", throw_error=true)
    err = ccall(func_ptr,
                Int32,
                (Float64, Float64, Float64, Float64, Float64, Float64, Float64),
                R, E, B, D, G, H, Q);

end

main()
