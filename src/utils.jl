using LinearAlgebra

function skew(v)
    [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end

function matrix_isapprox(A,B;atol = 1e-8)
    # @test size(A) == size(B)
    # @test isapprox(vec(A),vec(B),atol = atol)
    if (size(A) == size(B)) && isapprox(vec(A),vec(B),atol = atol)
        return true
    else
        return false
    end
end
