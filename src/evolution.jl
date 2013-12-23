export ## Methods
       unitary_propagator


const liblapack = Base.liblapack_name
import Base.LinAlg: BlasChar, BlasInt, blas_int

function expm_eigen(A::Matrix, t)
    #Calculates exp(t*A) via eigenvalue decomposition and assuming Hermitian matrix
    F = eigfact(Hermitian(A))

    # V * diagm(exp(t*D)) * V'
    return scale(F[:vectors], exp(t*F[:values])) * F[:vectors]'
end

function allocate_workspace(dim::Int)
    #Allocate workspace for Hermitian eigensolver for matrices of dimension dim

    #Create a random hermitian matrix
    A = rand((dim,dim)) + 1im*rand((dim,dim))
    H = A + A'


    #Make the initial call as a workspace query with lwork=-1 
        # subroutine zheevr   ( 
        #       character   JOBZ,
        #       character   RANGE,
        #       character   UPLO,
        #       integer     N,
        #       complex*16, dimension( lda, * )     A,
        #       integer     LDA,
        #       double precision    VL,
        #       double precision    VU,
        #       integer     IL,
        #       integer     IU,
        #       double precision    ABSTOL,
        #       integer     M,
        #       double precision, dimension( * )    W,
        #       complex*16, dimension( ldz, * )     Z,
        #       integer     LDZ,
        #       integer, dimension( * )     ISUPPZ,
        #       complex*16, dimension( * )  WORK,
        #       integer     LWORK,
        #       double precision, dimension( * )    RWORK,
        #       integer     LRWORK,
        #       integer, dimension( * )     IWORK,
        #       integer     LIWORK,
        #       integer     INFO )   

    jobz = 'V' # compute both eigenvalues and eigenvectors
    range = 'A' # all eigenvalues will be found
    uplo = 'U' # upper triangle storage (but we have full matrices anyways)
    n = blas_int(dim)
    vl = 0.0 #lower limit of eigenvalues: range = 'A'
    vu = 0.0 #upper limit of eigenvalues: range = 'A'
    il = blas_int(0) #lower number of eigenvalue to find: range = 'A'
    iu = blas_int(0) #upper number of eigenvalue to find: range = 'A'
    abstol = -1.0 # precision: negative for full precision
    m = Array(BlasInt, 1) #number of non-zero eigenvalues found 
    w = zeros(Float64, dim) #array for eigenvalues
    z = zeros(Complex128, (dim,dim)) #array for eigenvectors
    ldz = blas_int(dim) # number of rows in Z

    #Make initial allocation of empty work arrays
    isuppz = Array(BlasInt, 2*dim)
    work  = Array(Complex128, 1)
    lwork = blas_int(-1)
    rwork = Array(Float64, 1)
    lrwork = blas_int(-1)
    iwork = Array(BlasInt, 1)
    liwork = blas_int(-1)
    info  = Array(BlasInt, 1)

    ccall(("zheevr_", liblapack), Void,
                    (Ptr{BlasChar}, Ptr{BlasChar}, Ptr{BlasChar}, Ptr{BlasInt},
                     Ptr{Complex128}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
                     Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Complex128}, Ptr{BlasInt},
                     Ptr{Float64}, Ptr{Complex128}, Ptr{BlasInt}, Ptr{BlasInt},
                     Ptr{Complex128}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                     Ptr{BlasInt},Ptr{BlasInt}, Ptr{BlasInt}),
                    &jobz, &range, &uplo, &n, 
                    A, &n, &vl, &vu, 
                    &il, &iu, &abstol, m,
                    w, z, &ldz, isuppz,
                    work, &lwork, rwork, &lrwork,
                    iwork, &liwork, info)

    #Update the work arrays 
    lwork = blas_int(real(work[1]))
    work = Array(Complex128, lwork)
    lrwork = blas_int(rwork[1])
    rwork = Array(Float64, lrwork)
    liwork = iwork[1]
    iwork = Array(BlasInt, liwork)


    #Return all the allocated memory
    return (jobz, range, uplo, n, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)

end

function expm_eigen!(A::Matrix, t, jobz, range, uplo, n, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
    #Matrix exponential exp(t*A) for Hermitian matrix
    #Works via Hermitian eigensolver
    #This version requires preallocated workspace memory and destroys the input matrix A

    #Directly call the LAPACK function 
    ccall(("zheevr_", liblapack), Void,
                    (Ptr{BlasChar}, Ptr{BlasChar}, Ptr{BlasChar}, Ptr{BlasInt},
                     Ptr{Complex128}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
                     Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Complex128}, Ptr{BlasInt},
                     Ptr{Float64}, Ptr{Complex128}, Ptr{BlasInt}, Ptr{BlasInt},
                     Ptr{Complex128}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                     Ptr{BlasInt},Ptr{BlasInt}, Ptr{BlasInt}),
                    &jobz, &range, &uplo, &n, 
                    A, &n, &vl, &vu, 
                    &il, &iu, &abstol, m,
                    w, z, &ldz, isuppz,
                    work, &lwork, rwork, &lrwork,
                    iwork, &liwork, info)

    #TODO: should we be checking m to make sure there are dim eigenvalues?
    return scale(z, exp(t*w)) * z'

end

function unitary_propagator(sys::CompositeQSystem, timeStep::Float64, startTime::Float64, endTime::Float64; parallelize=true)

    #Preallocate Hamiltonian memory
    Ham = zeros(Complex128, (dim(sys), dim(sys)))

    #Allocate workspace for the matrix exponential
    workspace = allocate_workspace(dim(sys))

    if parallelize
        Uprop = @parallel (*) for time = startTime:timeStep:endTime
            #a *= b expands to a = a*b
            hamiltonian_add!(Ham, sys, time)
            # expm_eigen(Ham, 1im*2pi*timeStep)
            expm_eigen!(Ham, 1im*2pi*timeStep, workspace...)
        end
    else
        Uprop = eye(dim(sys))
        for time = startTime:timeStep:endTime
            #a *= b expands to a = a*b
            hamiltonian_add!(Ham, sys, time)
            # expm_eigen(Ham, 1im*2pi*timeStep)
            Uprop *= expm_eigen!(Ham, 1im*2pi*timeStep, workspace...)
        end
    end
    return Uprop'
end


