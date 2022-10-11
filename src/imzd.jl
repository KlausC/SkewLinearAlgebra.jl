
function imzd(v::AbstractVector, SKEW::Bool, Z)
    E = [0; v]
    ierr = imzd!(E, SKEW, Z)
    ierr != 0 && error("did not converge at $ierr")
    E
end

"""
    imzd!(S::vector, skew::Bool[, Z::matrix])

Return the eigenvalues of a symmetric or skew symmetric real tridiagonal
matrix with zero diagonal. The subdiagonal is given by `S`.
The optional matrix `Z` of appropriate size is overwritten with
the eigenvectors.

### Arguments
- `S`: real vector defining the subdiagonal
- `skew`: if skew symmetric, else symmetric
- `Z`: `nothing` if no eigenvectors, else of type `AbstractMatrix`

# Extended help

The original documentation of FORTRAN77 source code.

```
SUBROUTINE IMZD ( N, E, MATZ, SKEW, NZ, Z, IERR )                 IMZD0001
C
C ****
C
C  FUNCTION  -  COMPUTE THE EIGENVALUES AND OPTIONALLY THE EIGENVECTORS
C                 OF A SYMMETRIC TRIDIAGONAL MATRIX WITH ZERO DIAGONALS
C                 OR A SKEW-SYMMETRIC TRIDIAGONAL MATRIX USING AN
C                 IMPLICIT QR-TYPE ITERATION
C
C  PARAMETERS
C
C     N        - INPUT INTEGER SPECIFYING THE ORDER OF THE TRIDIAGONAL
C                  MATRIX
C
C     E(N)     - ON INPUT, ARRAY CONTAINING THE LOWER SUBDIAGONAL
C                  ELEMENTS OF THE TRIDIAGONAL MATRIX IN ITS LAST N-1
C                  POSITIONS. E(1) IS ARBITRARY.
C                ON OUTPUT, ARRAY CONTAINS THE EIGENVALUES. THE NON-ZERO
C                  EIGENVALUES OCCUR IN PAIRS WITH OPPOSITE SIGNS AND
C                  ARE FOUND IN ADJACENT LOCATIONS IN E. THE EIGENVALUES
C                  OF SYMMETRIC MATRICES ARE REAL AND THE EIGENVALUES
C                  OF SKEW-SYMMETRIC MATRICES ARE PURELY IMAGINARY
C                  COMPLEX NUMBERS. IF AN ERROR EXIT IS MADE, THE
C                  EIGENVALUES ARE CORRECT FOR INDICES IERR+1,IERR+2...N
C
C     MATZ     - INPUT LOGICAL VARIABLE SPECIFYING THE EIGENVECTOR
C                  OPTION
C                  = .TRUE.   EIGENVECTORS ARE TO BE COMPUTED
C                  = .FALSE.  EIGENVECTORS ARE NOT TO BE COMPUTED
C
C     SKEW     - INPUT LOGICAL VARIABLE SPECIFYING TYPE OF INPUT MATRIX
C                  = .TRUE.   INPUT TRIDIAGONAL MATRIX IS SKEW-SYMMETRIC
C                  = .FALSE.  INPUT TRIDIAGONAL MATRIX IS SYMMETRIC WITH
C                               ZERO DIAGONALS
C                  SKEW IS NOT REFERENCED IF MATZ = .FALSE.
C
C     NZ       - INPUT INTEGER SPECIFYING THE ROW DIMENSION OF Z AS
C                  DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT
C
C     Z(NZ,N)  - OUTPUT ARRAY CONTAINING THE ORTHOGONAL EIGENVECTORS
C                  OF THE INPUT TRIDIAGONAL MATRIX. EIGENVECTORS CORRE-
C                  SPONDING TO ZERO EIGENVALUES ARE NORMALIZE TO UNIT
C                  2-NORM (LENGTH) AND THOSE CORRESPONDING TO NON-ZERO
C                  EIGENVALUES HAVE 2-NORM OF SQUARE ROOT 2. IF THE J-TH
C                  EIGENVALUE IS ZERO OR REAL (I.E. E(J)), ITS EIGEN-
C                  VECTOR IS FOUND IN THE J-TH COLUMN OF Z. IF THE J-TH
C                  EIGENVALUE IS IMAGINARY (I.E. E(J)*I) WITH E(J+1) =
C                  -E(J), THE REAL PART OF ITS EIGENVECTOR IS FOUND IN
C                  THE J-TH COLUMN OF Z AND ITS IMAGINARY PART FOUND IN
C                  THE (J+1)-TH COLUMN. IF AN ERROR EXIT IS MADE, Z
C                  CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C                  EIGENVALUES.
C                  Z IS NOT REFERENCED IF MATZ = .FALSE.
C
C     IERR     - OUTPUT ERROR CODE
C                  = 0   NORMAL RETURN (ALL EIGENVALUES/VECTORS FOUND)
C                  = J   IF THE J-TH EIGENVALUE HAS NOT BEEN DETERMINED
C                          AFTER 30 ITERATIONS
C
C  REQUIRED FUNCTIONS - ABS,SIGN,SQRT,MOD
C
C ****
```
"""
function imzd!(E::AbstractVector{T}, SKEW::Bool, Z::Union{Nothing,AbstractMatrix{T}}) where {T<:Real}
    IERR = 0
    N = length(E)
    MATZ = Z !== nothing
    ZERO = zero(T)
    ONE = one(T)

    if MATZ
        LinearAlgebra.checksquare(Z)
        size(Z, 1) == N || throw(DimensionMismatch("matrix vs. eigenvectors"))
        # *** PLACE IDENTITY MATRIX IN Z
        for I = 1:N, J = 1:N
            Z[I, J] = ifelse(I == J, ONE, ZERO)
        end
    end
    M = N
    E[1] = ZERO
    ITS = 0

    while M >= 2 && IERR == 0
        M0 = M
        # *** SEARCH FOR NEXT SUBMATRIX TO SOLVE  (MATRIX SPLITTING)
        F = ZERO
        L0 = 2
        for J = M-1:-1:1
            G = abs(E[J+1])
            TMAG = abs(E[J]) + F
            TEST = TMAG + G
            if TEST == TMAG
                L0 += J
                break
            end
            F = G
        end
        L0 - 1 == M && @goto CONVERGED
        if MATZ && SKEW
            # *** PLACE CORRECT SIGN ON IDENTITY DIAGONALS
            for I = L0-1:4:M
                Z[I, I] = -Z[I, I]
                if I+3 <= M
                    Z[I+3, I+3] = -Z[I+3, I+3]
                end
            end
        end
        L0 == M && @goto CONVERGED
        L = L0
        iseven(M - L0) && @goto L230
        # *** FIND ZERO EIGENVALUE OF ODD ORDERED SUBMATRICES
        C = ZERO
        S = -ONE
        for K = M-1:-2:L0
            Q = -S * E[K+1]
            E[K+1] = C * E[K+1]
            if abs(E[K]) <= abs(Q)
                C = E[K] / Q
                R = sqrt(C * C + ONE)
                E[K] = Q * R
                S = ONE / R
                C = C * S
            else
                S = Q / E[K]
                R = sqrt(ONE + S * S)
                E[K] = E[K] * R
                C = ONE / R
                S = S * C
            end
            if MATZ
                # *** ACCUMULATE TRANSFORMATIONS FOR EIGENVECTORS
                Z[K-1, M] = -S * Z[K-1, K-1]
                Z[K-1, K-1] = C * Z[K-1, K-1]
                for J = K+1:2:M
                    Z[J, K-1] = S * Z[J, M]
                    Z[J, M] = C * Z[J, M]
                end
            end
        end
        M = M - 1
        L0 == M && @goto CONVERGED
        @label L200
        while true
            # *** CHECK FOR CONVERGENCE OR SMALL SUBDIAGONAL ELEMENT
            L = L0
            for K = M-1:-2:L0
                TMAG = abs(E[K+1]) + abs(E[K-1])
                TEST = TMAG + E[K]
                if TEST == TMAG
                    L = K + 1
                    break
                end
            end
            L == M && @goto CONVERGED

            @label L230
            ITS = ITS + 1
            if ITS > 300
                # *** ERROR EXIT
                IERR = M
                break
            end
            # *** FORM SHIFT
            F = E[M-3]
            G = E[M-2]
            C = E[M-1]
            S = E[M]
            P = ((C - F) * (C + F) + (S - G) * (S + G)) / (2.0 * G * C)
            R = sqrt(P * P + ONE)
            Q = (G / (P + copysign(R, P))) - C
            F = E[L]
            E[L-1] = ((F - S) * (F + S) + C * Q) / F
            # *** PERFORM ONE IMPLICIT QR ITERATION ON CHOLESKY FACTOR
            LS = L0 - 1
            C = ONE
            S = ONE
            for I = L:M-1
                Q = S * E[I+1]
                E[I+1] = C * E[I+1]
                if abs(E[I-1]) <= abs(Q)
                    C = E[I-1] / Q
                    R = sqrt(C * C + ONE)
                    E[I-1] = Q * R
                    S = ONE / R
                    C = C * S
                else
                    S = Q / E[I-1]
                    R = sqrt(ONE + S * S)
                    E[I-1] = E[I-1] * R
                    C = ONE / R
                    S = S * C
                end
                F = E[I+1]
                E[I+1] = -S * E[I] + C * F
                E[I] = C * E[I] + S * F
                # *** ACCUMULATE TRANSFORMATIONS FOR EIGENVECTORS
                if MATZ
                    for J = LS:2:M0
                        F = Z[J, I+1]
                        Z[J, I+1] = -S * Z[J, I-1] + C * F
                        Z[J, I-1] = C * Z[J, I-1] + S * F
                    end
                    LS = LS == L0 - 1 ? L0 : L0 - 1
                end
            end
            E[L-1] = ZERO
        end
        @label CONVERGED
        ## println("converged $ITS $L0:$M, $(E[max(1, L0-2, M-4):M])")
        ITS = 0
        if L0 - 1 == M
            # *** ITERATION CONVERGED TO ONE ZERO EIGENVALUE
            E[M] = ZERO
            M = M - 1
        else
            # *** ITERATION CONVERGED TO EIGENVALUE PAIR
            while true
                E[M-1] = E[M]
                E[M] = -E[M]
                M = M - 2
                M != L0 && break
            end
        end

        M > L0 && @goto L200
        if MATZ && !SKEW
            # *** COMPUTE EIGENVECTORS FROM ORTHONORMAL COLUMNS OF Z IF NOT SKEW
            K = M0
            while K >= L0 - 1
                if E[K] != ZERO
                    for J = L0-1:2:M0
                        Z[J, K] = Z[J, K-1]
                        F = Z[J+1, K]
                        Z[J+1, K-1] = F
                        Z[J+1, K] = -F
                    end
                    K = K-1
                end
                K = K - 1
            end
        end
    end
    return IERR
end
