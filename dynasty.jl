
using Libdl
import Printf.@printf

function main()
    R = 0.07
    E = 2.0
    B = 0.17
    D = 0.42
    G = 0.09
    H = 0.1
    Q = 0.4

    pars    = [R, E, B, D, G, H, Q]
    y0      = [0.1, 0.1, 0.8];
    ttran   = 10e3
    tend    = 50e3
    max_nev = 10
    atol    = [1e-10, 1e-10, 1e-10]
    rtol    = Ref{Float64}(1e-12)
    sol     = zeros(max_nev * 4)
    
    lib = dlopen("dynasty.dylib", RTLD_LAZY | RTLD_GLOBAL)
    func_ptr = dlsym(lib, "integrate", throw_error=true)
    @printf("C:\n")
    nev = ccall(func_ptr,
                Int32,
                (Ptr{Float64}, Ptr{Float64}, Float64, Float64,
                 Int, Ptr{Float64}, Ref{Float64}, Ptr{Float64}),
                pars, y0, ttran, tend, max_nev, atol, rtol, sol);

    sol = reshape(sol, 4, :)
    @printf("\nJulia:\n")
    for i in 1 : nev
        @printf("[%03zu/%03zu] %12.4f %13.10f %13.10f %13.10f\n",
                i, nev, sol[1,i], sol[2,i], sol[3,i], sol[4,i])
    end
end

main()
