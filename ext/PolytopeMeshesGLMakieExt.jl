module PolytopeMeshesGLMakieExt

using PolytopeMeshes
using GLMakie

PolytopeMeshes.backend_is_interactive(::Val{:GLMakie}) = true

# Open each interactive viewer in its own window
function PolytopeMeshes.display_interactive(fig, ::Val{:GLMakie})
    if get(ENV, "CI", "false") == "false"
        display(GLMakie.Screen(), fig)
    end
    return fig
end

end # module
