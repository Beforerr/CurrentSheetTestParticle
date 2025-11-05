using FHist

const DEFAULT_W_BINEDGE = range(-1, 1, length = 46)
const DEFAULT_W_BINEDGES = (DEFAULT_W_BINEDGE, DEFAULT_W_BINEDGE)

"""
right stochastic matrix, with each row summing to 1., m〚i,j〛 = Probability[x[k+1]=j|x[k]=i]
"""
function transition_matrix!(h::Hist2D)
    m = bincounts(h)
    return m ./= sum(m, dims = 2)
end

function transition_matrix(df, xsym = :μ0, ysym = :μ1; binedges = nothing, kw...)
    x = df[!, xsym]
    y = df[!, ysym]
    binedges = @something binedges begin
        xedges = sort(unique(x))
        (xedges, xedges)
    end
    h = Hist2D((x, y); binedges, kw...)
    mat = transition_matrix!(h)
    bcenters = bincenters(h)
    dims = (Dim{xsym}(bcenters[1]), Dim{ysym}(bcenters[2]))
    return DimArray(mat, dims)
end
