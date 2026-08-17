"""Custom, minimal depth-diffusion Step for the N-layer biofilm model (see
composites/eps_biofilm_depth_scaffold.py for the full design rationale). Built
as a small local replacement for spatio_flux's DiffusionAdvection after a direct smoke test
confirmed a real shape-order bug for any non-square `n_bins` grid: its
inputs()/outputs() declare `_shape` from `n_bins` in (nx,ny) order, but its
internal update() enforces the array shape in (ny,nx) order -- these only ever
agree when nx==ny (true of every existing spatio_flux test/composite, none of
which use a thin grid). `external/spatio-flux` is a pinned git submodule, not
vendored/owned code -- rather than patch it, this file reimplements the SAME
physics (explicit finite-difference diffusion with ghost-cell dirichlet/
neumann boundaries, matching DiffusionAdvection's own convention) directly for
the simpler 1D case actually needed here: N depth layers, one dirichlet
boundary (bulk fluid, at the surface/top layer) and one neumann boundary
(glass surface, at the basal/bottom layer). No advection term -- horizontal
chamber flow doesn't move nutrients vertically through the matrix (see
composites/eps_biofilm_depth_scaffold.py for the full rationale).
"""
import math

from process_bigraph import Process


class DepthDiffusionStep(Process):
    """Diffuses each amino acid's concentration across `num_layers` depth
    layers. State shape: `substrates[aa_id][str(layer_index)] = concentration
    (mM)`. Layer index 0 = surface (dirichlet at `bulk_concentration[aa_id]`),
    layer index `num_layers - 1` = basal (neumann / no-flux) -- this
    convention is internal to this Step, chosen to match
    `eps_biofilm_depth_scaffold.py`'s own layer-index convention directly, no
    row-order translation needed anywhere.
    """

    config_schema = {
        "num_layers": {"_type": "integer", "_default": 3},
        "aa_ids": "list[string]",
        "diffusion_coeff_cm2_per_h": {"_type": "float", "_default": 0.009},
        "depth_per_layer_cm": "float",
        "bulk_concentration": "map[float]",
        # Explicit-Euler 1D diffusion stability limit is dt <= dz^2/(2D); this
        # Step sub-steps internally (like DiffusionAdvection's own
        # _compute_stable_dt) so the outer composite's dt never has to respect
        # this by itself. cfl_safety<1 leaves headroom below the hard limit.
        "cfl_safety": {"_type": "float", "_default": 0.4},
    }

    def initialize(self, config):
        pass

    def inputs(self):
        return {
            "substrates": {
                "_type": "map",
                "_value": {"_type": "map", "_value": {"_type": "concentration", "_units": "mM"}},
            },
        }

    def outputs(self):
        return {
            "substrates": {
                "_type": "map",
                "_value": {"_type": "map", "_value": "concentration_delta"},
            },
        }

    def update(self, inputs, interval):
        n = self.config["num_layers"]
        dz = self.config["depth_per_layer_cm"]
        D = self.config["diffusion_coeff_cm2_per_h"]
        bulk = self.config["bulk_concentration"]
        cfl_safety = self.config["cfl_safety"]

        dt_stable = cfl_safety * dz * dz / (2.0 * D)
        n_steps = max(1, math.ceil(interval / dt_stable))
        sub_dt = interval / n_steps

        substrates_in = inputs["substrates"]
        deltas = {aa_id: {str(i): 0.0 for i in range(n)} for aa_id in substrates_in}

        for aa_id, layer_map in substrates_in.items():
            c = [float(layer_map[str(i)]) for i in range(n)]
            c0 = list(c)
            top_val = bulk.get(aa_id, c[0])

            for _ in range(n_steps):
                # Ghost-padded array (length n+2): g[0] = bottom ghost (neumann,
                # mirrors the adjacent interior cell -- zero flux at the glass
                # surface); g[n+1] = top ghost (dirichlet, reflected around
                # top_val -- matches DiffusionAdvection's own g = 2*val - c
                # convention, so this is a real, established boundary formula,
                # not an ad hoc one). g[1:n+1] holds the interior layers.
                g = [0.0] * (n + 2)
                g[1:n + 1] = c
                g[0] = c[0]
                g[n + 1] = 2.0 * top_val - c[n - 1]

                new_c = list(c)
                for i in range(n):
                    laplacian = (g[i] - 2.0 * g[i + 1] + g[i + 2]) / (dz * dz)
                    new_c[i] = max(0.0, c[i] + sub_dt * D * laplacian)
                c = new_c

            for i in range(n):
                deltas[aa_id][str(i)] = c[i] - c0[i]

        return {"substrates": deltas}
