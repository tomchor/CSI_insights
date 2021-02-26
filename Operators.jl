#module Operators

using Oceananigans.Operators
using KernelAbstractions: @index, @kernel
using Oceananigans.Fields: location
using Oceananigans.AbstractOperations: AbstractOperation


#++++ Simple derivates
@kernel function ∂xᶜ!(result, grid, operand)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ∂xᶜᵃᵃ(i, j, k, grid, operand) # F, A, A  → C, A, A
end

@kernel function ∂xᶠ!(result, grid, operand)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ∂xᶠᵃᵃ(i, j, k, grid, operand) # C, A, A  → F, A, A
end



@kernel function ∂yᶜ!(result, grid, operand)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ∂yᵃᶜᵃ(i, j, k, grid, operand) # A, F, A  → A, C, A
end

@kernel function ∂yᶠ!(result, grid, operand)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ∂yᵃᶠᵃ(i, j, k, grid, operand) # A, C, A  → A, F, A
end



@kernel function ∂zᶜ!(result, grid, operand)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ∂zᵃᵃᶜ(i, j, k, grid, operand) # A, A, F  → A, A, C
end

@kernel function ∂zᶠ!(result, grid, operand)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ∂zᵃᵃᶠ(i, j, k, grid, operand) # A, A, C  → A, A, F
end
#-----



#++++ Multiplication times a simple derivate
@kernel function ϕ∂xᶜ!(result, grid, operand, ϕ)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ϕ * ∂xᶜᵃᵃ(i, j, k, grid, operand) # F, ϕ, A  → C, A, A
end

@kernel function ϕ∂xᶠ!(result, grid, operand, ϕ)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ϕ * ∂xᶠᵃᵃ(i, j, k, grid, operand) # C, A, A  → F, A, A
end



@kernel function ϕ∂yᶜ!(result, grid, operand, ϕ)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ϕ * ∂yᵃᶜᵃ(i, j, k, grid, operand) # A, F, A  → A, C, A
end

@kernel function ϕ∂yᶠ!(result, grid, operand, ϕ)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ϕ * ∂yᵃᶠᵃ(i, j, k, grid, operand) # A, C, A  → A, F, A
end



@kernel function ϕ∂zᶜ!(result, grid, operand, ϕ)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ϕ * ∂zᵃᵃᶜ(i, j, k, grid, operand) # A, A, F  → A, A, C
end

@kernel function ϕ∂zᶠ!(result, grid, operand, ϕ)
    i, j, k = @index(Global, NTuple)
    @inbounds result[i, j, k] = ϕ * ∂zᵃᵃᶠ(i, j, k, grid, operand) # A, A, C  → A, A, F
end
#-----



#end # module

