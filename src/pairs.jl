@inline pair_dimension(rank::Int) = rank * (rank + 1) ÷ 2
@inline pair_index(first::Int, second::Int) = begin
    left, right = minmax(first, second)
    right * (right - 1) ÷ 2 + left
end
@inline pair_scale(first::Int, second::Int) =
    first == second ? 1.0 : sqrt(2.0)

function _pair_transform!(output, orbital_map, old_rank::Int, new_rank::Int)
    old_pairs = pair_dimension(old_rank)
    new_pairs = pair_dimension(new_rank)
    size(output, 1) >= old_pairs && size(output, 2) >= new_pairs ||
        throw(DimensionMismatch("pair-transform capacity is too small"))
    fill!(@view(output[1:old_pairs, 1:new_pairs]), 0.0)
    @inbounds for second_new = 1:new_rank, first_new = 1:second_new
        new_pair = pair_index(first_new, second_new)
        new_scale = pair_scale(first_new, second_new)
        for second_old = 1:old_rank, first_old = 1:second_old
            old_pair = pair_index(first_old, second_old)
            value = orbital_map[first_old, first_new] *
                orbital_map[second_old, second_new]
            first_old != second_old && (value +=
                orbital_map[second_old, first_new] *
                orbital_map[first_old, second_new])
            output[old_pair, new_pair] = new_scale * value /
                pair_scale(first_old, second_old)
        end
    end
    @view output[1:old_pairs, 1:new_pairs]
end

@inline function _raw_integral(pair_field, a, c, b, d)
    pair_field[pair_index(a, c), pair_index(b, d)] /
        (pair_scale(a, c) * pair_scale(b, d))
end
