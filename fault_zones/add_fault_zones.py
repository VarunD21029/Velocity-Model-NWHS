#!/usr/bin/env python3
"""
Imprint fault damage zones on the SPECFEM3D tomography model of the NW Himalaya.

For every fault in the fault file (surface trace + dip profile), each grid node of
tomography_model.xyz gets a Vs reduction that is

    * strongest on the fault (default 10 % at the surface),
    * carried down-dip along the fault plane, following a ramp-flat / listric
      geometry, with a depth-dependent amplitude (default 10 % -> 5 % at 20 km),
    * tapered smoothly to zero away from the fault plane and beyond the trace ends.

The Main Himalayan Thrust (MHT) is modelled as one flat-ramp-flat decollement
(MFT frontal ramp -> upper flat -> mid-crustal ramp -> lower flat), and the splay
thrusts (JMT, MBT, MCT) root into it, so the thrusts form one connected system.
Where zones overlap, the largest reduction is used, so reductions never stack.

The grid order, coordinates and number formatting of the input file are kept
unchanged; only the Vs column (and optionally Vp / rho) and the min/max values in
the header are modified.

Usage:
    python3 add_fault_zones.py \
        --input ../velocity_model_file.zip \
        --faults faults_nwh.json \
        --output tomography_model_faults.xyz \
        [--figures figures/]

Only numpy is needed; matplotlib is optional and only used for --figures.
"""

import argparse
import io
import json
import math
import os
import zipfile

import numpy as np


# --------------------------------------------------------------------------
# WGS84 -> UTM (Krueger series, accurate to well below 1 m within a zone)
# --------------------------------------------------------------------------

def lonlat_to_utm(lon, lat, zone):
    a = 6378137.0
    f = 1.0 / 298.257223563
    k0 = 0.9996
    n = f / (2.0 - f)
    A = a / (1.0 + n) * (1.0 + n**2 / 4.0 + n**4 / 64.0)
    alpha = (n / 2.0 - 2.0 * n**2 / 3.0 + 5.0 * n**3 / 16.0,
             13.0 * n**2 / 48.0 - 3.0 * n**3 / 5.0,
             61.0 * n**3 / 240.0)
    lon = np.radians(np.asarray(lon, dtype=float))
    lat = np.radians(np.asarray(lat, dtype=float))
    lon0 = np.radians(zone * 6.0 - 183.0)
    e = math.sqrt(f * (2.0 - f))
    t = np.sinh(np.arctanh(np.sin(lat)) - e * np.arctanh(e * np.sin(lat)))
    xi_p = np.arctan2(t, np.cos(lon - lon0))
    eta_p = np.arctanh(np.sin(lon - lon0) / np.sqrt(1.0 + t**2))
    xi, eta = xi_p.copy(), eta_p.copy()
    for j, aj in enumerate(alpha, start=1):
        xi += aj * np.sin(2 * j * xi_p) * np.cosh(2 * j * eta_p)
        eta += aj * np.cos(2 * j * xi_p) * np.sinh(2 * j * eta_p)
    easting = 500000.0 + k0 * A * eta
    northing = k0 * A * xi          # northern hemisphere
    return easting, northing


# --------------------------------------------------------------------------
# Tomography file I/O
# --------------------------------------------------------------------------

def read_tomography(path):
    """Return (header_lines, data_lines, array[N, 6]) for a .xyz or a .zip containing one."""
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            name = [n for n in zf.namelist() if n.endswith(".xyz")][0]
            text = zf.read(name).decode()
    else:
        with open(path) as fh:
            text = fh.read()
    lines = text.splitlines()
    # SPECFEM header: comment lines + 4 numeric lines; data starts after the last
    # comment line that precedes the first record.
    first_data = None
    numeric_seen = 0
    for i, line in enumerate(lines):
        if line.lstrip().startswith("#") or not line.strip():
            continue
        numeric_seen += 1
        if numeric_seen == 5:
            first_data = i
            break
    header = lines[:first_data]
    while first_data < len(lines) and lines[first_data].lstrip().startswith("#"):
        header.append(lines[first_data])
        first_data += 1
    data = np.loadtxt(io.StringIO("\n".join(lines[first_data:])))
    return header, data


def write_tomography(path, header, data):
    vp, vs, rho = data[:, 3], data[:, 4], data[:, 5]
    out_header = list(header)
    # The 4th numeric header line holds vp/vs/rho min & max.
    numeric_idx = [i for i, l in enumerate(out_header) if l.strip() and not l.lstrip().startswith("#")]
    out_header[numeric_idx[3]] = "%.2f %.2f %.2f %.2f %.2f %.2f " % (
        vp.min(), vp.max(), vs.min(), vs.max(), rho.min(), rho.max())
    with open(path, "w") as fh:
        fh.write("\n".join(out_header) + "\n")
        np.savetxt(fh, data, fmt="%.2f", delimiter=" ", newline=" \n")


# --------------------------------------------------------------------------
# Fault geometry
# --------------------------------------------------------------------------

class Fault:
    def __init__(self, spec, zone):
        self.name = spec["name"]
        self.short = spec.get("short", self.name)
        self.mht = spec.get("mht_geometry")          # only for the Main Himalayan Thrust
        self.is_mht = self.mht is not None
        self.terminate_at_mht = spec.get("terminate_at_mht", False)
        self.max_depth = spec["max_depth_km"] * 1e3
        self.trace_lonlat = np.array(spec["trace"], dtype=float)
        x, y = lonlat_to_utm(self.trace_lonlat[:, 0], self.trace_lonlat[:, 1], zone)
        tr = np.column_stack([x, y])

        # Orient trace so the fault dips to the right of the strike direction
        # (Aki & Richards convention); the right-hand normal then points down-dip.
        seg = np.diff(tr, axis=0)
        right = np.column_stack([seg[:, 1], -seg[:, 0]])
        az = math.radians(spec["dip_azimuth_deg"])
        if np.sum(right @ np.array([math.sin(az), math.cos(az)])) < 0:
            tr = tr[::-1]
            self.trace_lonlat = self.trace_lonlat[::-1]
        self.trace = tr

        if not self.is_mht:
            # Fixed cross-section polyline (down-dip distance, depth) in metres. The
            # last dip continues well past max_depth; the fault is cut at max_depth
            # (or where it meets the MHT) later.
            prof = spec["dip_profile"]
            d, z = [0.0], [0.0]
            for k, (ztop, dip) in enumerate(prof):
                zbot = prof[k + 1][0] if k + 1 < len(prof) else max(self.max_depth / 1e3, ztop) + 60.0
                d.append(d[-1] + (zbot - ztop) * 1e3 / math.tan(math.radians(dip)))
                z.append(zbot * 1e3)
            self.sec_d = np.array(d)
            self.sec_z = np.array(z)
            self.sec_dip = np.radians([p[1] for p in prof])

    # -- map view ---------------------------------------------------------
    def horizontal_coords(self, px, py):
        """Signed down-dip distance from the trace (+ = hanging wall side) and
        along-strike overshoot beyond the trace ends (0 inside)."""
        tr = self.trace
        best = np.full(px.shape, np.inf)
        dist_signed = np.zeros(px.shape)
        overshoot = np.zeros(px.shape)
        nseg = len(tr) - 1
        for i in range(nseg):
            a, b = tr[i], tr[i + 1]
            ab = b - a
            L2 = ab @ ab
            L = math.sqrt(L2)
            t_raw = ((px - a[0]) * ab[0] + (py - a[1]) * ab[1]) / L2
            t = np.clip(t_raw, 0.0, 1.0)
            cx, cy = a[0] + t * ab[0], a[1] + t * ab[1]
            dd = np.hypot(px - cx, py - cy)
            # signed perpendicular distance to the infinite line of this segment
            perp = ((px - a[0]) * ab[1] - (py - a[1]) * ab[0]) / L
            ov = np.zeros(px.shape)
            if i == 0:
                ov = np.maximum(ov, -t_raw * L)
            if i == nseg - 1:
                ov = np.maximum(ov, (t_raw - 1.0) * L)
            # beyond the ends use the perpendicular component (plane continues straight)
            dist_here = np.where(ov > 0, np.abs(perp), dd)
            closer = dist_here < best
            best = np.where(closer, dist_here, best)
            dist_signed = np.where(closer, np.sign(perp) * dist_here, dist_signed)
            overshoot = np.where(closer, ov, overshoot)
        return dist_signed, overshoot

    # -- cross-section ----------------------------------------------------
    def section(self, n, d_anchor=None, d_self=None):
        """Cross-section polyline for n map positions: vertex distances D and depths Z
        with shape (nvert, n), and segment dips with shape (nvert - 1, n).

        For the MHT the polyline is frontal ramp -> upper flat -> mid-crustal ramp ->
        lower flat. The top of the mid-crustal ramp is tied to the anchor fault (MCT):
        it lies ramp_top_offset_km down-dip of the anchor's surface trace, so the ramp
        follows the front of the Higher Himalaya rather than a fixed distance from the
        MFT (the Sub- and Lesser Himalaya are much wider in the Kangra re-entrant)."""
        if not self.is_mht:
            D = np.repeat(self.sec_d[:, None], n, axis=1)
            Z = np.repeat(self.sec_z[:, None], n, axis=1)
            DIP = np.repeat(self.sec_dip[:, None], n, axis=1)
            return D, Z, DIP

        g = self.mht
        rad = math.radians
        fr_dip, fr_bot = rad(g["frontal_ramp_dip_deg"]), g["frontal_ramp_bottom_km"] * 1e3
        uf_dip = rad(g["upper_flat_dip_deg"])
        rp_dip, rp_bot = rad(g["ramp_dip_deg"]), g["ramp_bottom_km"] * 1e3
        lf_dip = rad(g["lower_flat_dip_deg"])
        z_end = self.max_depth + 60e3

        d1 = fr_bot / math.tan(fr_dip)
        # MFT -> ramp-top distance measured along each position's own dip direction
        d_hinge = d_anchor - g["ramp_top_offset_km"] * 1e3
        Dh = np.maximum(d_self - d_hinge, d1 + 1e3)
        Zh = fr_bot + (Dh - d1) * math.tan(uf_dip)
        Zr = np.maximum(rp_bot, Zh + 1e3)
        Dr = Dh + (Zr - Zh) / math.tan(rp_dip)
        De = Dr + (z_end - Zr) / math.tan(lf_dip)
        one = np.ones(n)
        D = np.vstack([0 * one, d1 * one, Dh, Dr, De])
        Z = np.vstack([0 * one, fr_bot * one, Zh, Zr, z_end * one])
        DIP = np.vstack([fr_dip * one, uf_dip * one, rp_dip * one, lf_dip * one])
        return D, Z, DIP

    @staticmethod
    def depth_at(d, D, Z):
        """Depth of the fault plane at signed down-dip distance d (nan on the footwall
        side); D, Z are per-position polylines from section()."""
        zf = np.full(d.shape, np.nan)
        for i in range(D.shape[0] - 1):
            last = i == D.shape[0] - 2
            m = (d >= D[i]) & ((d < D[i + 1]) | last)
            slope = (Z[i + 1] - Z[i]) / (D[i + 1] - D[i])
            zf = np.where(m, Z[i] + (d - D[i]) * slope, zf)
        return np.where(d >= 0, zf, np.nan)

    @staticmethod
    def section_distance(d, z, D, Z, DIP):
        """Distance from points (d, z) to the fault polyline in the cross-section,
        the local dip at the closest point and the depth of that closest point."""
        best = np.full(d.shape, np.inf)
        dip = np.zeros(d.shape)
        zc_best = np.zeros(d.shape)
        for i in range(D.shape[0] - 1):
            ad, az_, bd, bz = D[i], Z[i], D[i + 1], Z[i + 1]
            abd, abz = bd - ad, bz - az_
            t = np.clip(((d - ad) * abd + (z - az_) * abz) / (abd**2 + abz**2), 0.0, 1.0)
            cd, cz = ad + t * abd, az_ + t * abz
            dist = np.hypot(d - cd, z - cz)
            closer = dist < best
            best = np.where(closer, dist, best)
            dip = np.where(closer, DIP[i], dip)
            zc_best = np.where(closer, cz, zc_best)
        return best, dip, zc_best


def cosine_taper(x, x0, x1):
    """1 for x <= x0, 0 for x >= x1, half-cosine in between."""
    w = np.clip((x - x0) / np.maximum(x1 - x0, 1e-9), 0.0, 1.0)
    return 0.5 * (1.0 + np.cos(np.pi * w))


# --------------------------------------------------------------------------
# Main computation
# --------------------------------------------------------------------------

def mht_setup(faults, uxy):
    """MHT fault and its cross-section polylines at the unique map positions."""
    mht = [f for f in faults if f.is_mht]
    if not mht:
        return None, None
    mht = mht[0]
    anchor = [f for f in faults if f.short == mht.mht["ramp_anchor"]][0]
    d_self, _ = mht.horizontal_coords(uxy[:, 0], uxy[:, 1])
    d_anchor, _ = anchor.horizontal_coords(uxy[:, 0], uxy[:, 1])
    return mht, mht.section(len(uxy), d_anchor=d_anchor, d_self=d_self) + (d_self,)


def compute_reduction(data, faults, settings, grid_spacing):
    x, y, zcoord = data[:, 0], data[:, 1], data[:, 2]
    depth = -zcoord                     # tomography z is elevation (m, <= 0)
    hx, hy, hz = grid_spacing
    h_horiz = max(hx, hy)

    W = settings["damage_zone_half_width_km"] * 1e3
    taper_factor = settings["taper_width_factor"]
    end_taper = settings["end_taper_km"] * 1e3
    margin = settings["splay_termination_margin_km"] * 1e3
    table = np.array(settings["vs_reduction_vs_depth"], dtype=float)
    amp_depth = np.interp(depth / 1e3, table[:, 0], table[:, 1]) / 100.0

    # Work on unique horizontal positions; every depth level repeats them.
    xy = np.column_stack([x, y])
    uxy, inv = np.unique(xy, axis=0, return_inverse=True)
    inv = inv.ravel()

    mht, mht_sec = mht_setup(faults, uxy)
    mht_depth = None
    if mht is not None:
        D, Z, _, d_u = mht_sec
        mht_depth = Fault.depth_at(d_u, D, Z)[inv]

    reduction = np.zeros(len(data))
    owner = np.full(len(data), -1)
    for k, f in enumerate(faults):
        d_u, ov_u = f.horizontal_coords(uxy[:, 0], uxy[:, 1])
        if f.is_mht:
            D_u, Z_u, DIP_u = mht_sec[:3]
        else:
            D_u, Z_u, DIP_u = f.section(len(uxy))
        d, ov = d_u[inv], ov_u[inv]
        D, Z, DIP = D_u[:, inv], Z_u[:, inv], DIP_u[:, inv]

        n, dip, zc = Fault.section_distance(d, depth, D, Z, DIP)
        # Make the zone at least as wide as the grid can resolve: half a horizontal
        # grid step across steep faults, half a vertical step across flat ones.
        w_eff = np.maximum.reduce([np.full_like(n, W),
                                   0.5 * h_horiz * np.sin(dip),
                                   0.5 * hz * np.cos(dip)])
        r = cosine_taper(n, w_eff, w_eff * (1.0 + taper_factor))
        r *= cosine_taper(ov, 0.0, end_taper)
        # down-dip limit of this fault (closest point below max depth is not on the fault)
        r *= cosine_taper(zc - f.max_depth, 0.0, W)

        if f.terminate_at_mht and mht_depth is not None:
            # splays root into the MHT: cut them where they reach it
            zf = Fault.depth_at(d, D, Z)
            beyond = np.nan_to_num(zf - mht_depth, nan=-np.inf)
            r *= cosine_taper(beyond, margin, margin + W)

        r *= amp_depth
        newer = r > reduction
        reduction = np.where(newer, r, reduction)
        owner = np.where(newer, k, owner)
    return reduction, owner


def apply(data, reduction, settings):
    out = data.copy()
    out[:, 4] *= 1.0 - reduction
    out[:, 3] *= 1.0 - settings["vp_reduction_ratio"] * reduction
    out[:, 5] *= 1.0 - settings["rho_reduction_ratio"] * reduction
    return out


def fault_surfaces(faults, xy, end_taper):
    """Depth (m) of every fault plane below the map positions xy (nan where the
    fault is absent), with splays cut where they meet the MHT. Also returns the
    MHT ramp-top depth at each position."""
    mht, mht_sec = mht_setup(faults, xy)
    surf, ramp_top = {}, None
    mht_depth = None
    if mht is not None:
        D, Z, _, d_m = mht_sec
        mht_depth = Fault.depth_at(d_m, D, Z)
        ramp_top = Z[2]
    for f in faults:
        d, ov = f.horizontal_coords(xy[:, 0], xy[:, 1])
        if f.is_mht:
            D, Z = mht_sec[:2]
        else:
            D, Z, _ = f.section(len(xy))
        zf = Fault.depth_at(d, D, Z)
        zf[(zf > f.max_depth) | (ov > end_taper)] = np.nan
        if f.terminate_at_mht and mht_depth is not None:
            zf[zf > mht_depth] = np.nan
        surf[f.short] = zf
    return surf, mht_depth, ramp_top


def summarize_geometry(faults, xy, end_taper):
    surf, mht_depth, ramp_top = fault_surfaces(faults, xy, end_taper)
    if mht_depth is None:
        return
    inside = np.isfinite(surf[[f.short for f in faults if f.is_mht][0]])
    print("MHT: ramp top at %.1f-%.1f km depth (10-90th percentile over the model)"
          % tuple(np.percentile(ramp_top[inside], [10, 90]) / 1e3))
    for f in faults:
        if not f.terminate_at_mht:
            continue
        zf = surf[f.short]
        junction = np.isfinite(zf) & (np.abs(zf - mht_depth) < 1.5e3)
        if junction.any():
            print("  %-4s roots into the MHT at %.1f-%.1f km depth"
                  % ((f.short,) + tuple(np.percentile(mht_depth[junction], [10, 90]) / 1e3)))
        else:
            print("  %-4s does NOT reach the MHT inside the model" % f.short)


# --------------------------------------------------------------------------
# Optional QC figures
# --------------------------------------------------------------------------

def make_figures(outdir, data, new, reduction, faults, header_dims, settings):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    nx, ny, nz = header_dims
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    pct = 100.0 * reduction
    colors = {"MFT": "tab:red", "JMT": "tab:orange", "MBT": "gold", "MCT": "tab:green",
              "STD": "tab:blue", "KCF": "tab:purple", "KKF": "tab:brown"}

    # Fault planes below one depth level of the grid (file order: y fastest, then x, then z)
    xy0 = data[:nx * ny, :2]
    surf, mht_depth, ramp_top = fault_surfaces(faults, xy0, settings["end_taper_km"] * 1e3)

    def label_traces(ax):
        for f in faults:
            ax.plot(f.trace[:, 0] / 1e3, f.trace[:, 1] / 1e3, "-", color=colors.get(f.short, "c"), lw=1.4)
            mid = f.trace[len(f.trace) // 2] / 1e3
            ax.text(mid[0], mid[1], f.short, color=colors.get(f.short, "c"), fontsize=8, weight="bold")
        ax.set_xlim(x.min() / 1e3, x.max() / 1e3)
        ax.set_ylim(y.min() / 1e3, y.max() / 1e3)
        ax.set_aspect("equal")
        ax.set_xlabel("UTM 43N easting (km)")
        ax.set_ylabel("northing (km)")

    # Map view of the reduction at a few depths
    depths = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0]
    zlev = np.unique(z)
    fig, axes = plt.subplots(2, 3, figsize=(18, 11), constrained_layout=True)
    for ax, dkm in zip(axes.ravel(), depths):
        zl = zlev[np.argmin(np.abs(zlev + dkm * 1e3))]
        m = z == zl
        sc = ax.scatter(x[m] / 1e3, y[m] / 1e3, c=pct[m], s=7, marker="s",
                        cmap="magma_r", vmin=0, vmax=10, linewidths=0)
        label_traces(ax)
        ax.set_title("Vs reduction at %.1f km depth" % (-zl / 1e3))
    fig.colorbar(sc, ax=axes, shrink=0.6, label="Vs reduction (%)")
    fig.savefig(os.path.join(outdir, "fault_zone_map_slices.png"), dpi=120)
    plt.close(fig)

    # MHT depth map (flat-ramp-flat geometry)
    if mht_depth is not None:
        X0 = xy0[:, 0].reshape(nx, ny) / 1e3
        Y0 = xy0[:, 1].reshape(nx, ny) / 1e3
        mfull = mht_depth.copy()
        mfull[~np.isfinite(surf["MFT"])] = np.nan
        fig, ax = plt.subplots(figsize=(9, 9), constrained_layout=True)
        cf = ax.contourf(X0, Y0, mfull.reshape(nx, ny) / 1e3, levels=np.arange(0, 32, 2), cmap="viridis_r")
        cs = ax.contour(X0, Y0, mfull.reshape(nx, ny) / 1e3, levels=np.arange(0, 32, 2), colors="k", linewidths=0.4)
        ax.clabel(cs, fmt="%d", fontsize=7)
        # surface projection of the mid-crustal ramp (between ramp top and ramp bottom)
        g = [f for f in faults if f.is_mht][0].mht
        on_ramp = (mfull >= ramp_top) & (mfull <= g["ramp_bottom_km"] * 1e3)
        ax.scatter(xy0[on_ramp, 0] / 1e3, xy0[on_ramp, 1] / 1e3, s=3, c="w", marker=".", label="mid-crustal ramp")
        label_traces(ax)
        ax.legend(loc="lower left")
        ax.set_title("Depth of the Main Himalayan Thrust (km)")
        fig.colorbar(cf, ax=ax, shrink=0.7, label="MHT depth (km)")
        fig.savefig(os.path.join(outdir, "mht_depth_map.png"), dpi=120)
        plt.close(fig)

    # Cross-sections with the fault planes overlaid
    def grid(col):
        return col.reshape(nz, nx, ny)

    X, Y, Z = grid(x), grid(y), grid(z)
    P, VS = grid(pct), grid(new[:, 4])
    fig, axes = plt.subplots(3, 2, figsize=(16, 12), constrained_layout=True)
    for row, ix in enumerate([nx // 5, nx // 2, 4 * nx // 5]):
        yy = Y[:, ix, :] / 1e3
        zz = -Z[:, ix, :] / 1e3
        a0 = axes[row, 0].pcolormesh(yy, zz, P[:, ix, :], cmap="magma_r", vmin=0, vmax=10, shading="auto")
        a1 = axes[row, 1].pcolormesh(yy, zz, VS[:, ix, :], cmap="viridis", shading="auto")
        line_y = Y[0, ix, :] / 1e3
        for ax in axes[row]:
            for f in faults:
                zf = surf[f.short].reshape(nx, ny)[ix] / 1e3
                ax.plot(line_y, zf, "-", color=colors.get(f.short, "c"), lw=1.0, label=f.short)
            ax.set_ylim(32, 0)
            ax.set_ylabel("depth (km)")
            ax.set_title("S-N section at easting %.0f km" % (X[0, ix, 0] / 1e3))
        fig.colorbar(a0, ax=axes[row, 0], label="Vs reduction (%)")
        fig.colorbar(a1, ax=axes[row, 1], label="Vs (m/s)")
    axes[0, 0].legend(loc="lower left", fontsize=7, ncol=4)
    for ax in axes[-1]:
        ax.set_xlabel("UTM northing (km)")
    fig.savefig(os.path.join(outdir, "fault_zone_cross_sections.png"), dpi=120)
    plt.close(fig)


# --------------------------------------------------------------------------

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", default=os.path.join(here, "..", "velocity_model_file.zip"))
    ap.add_argument("--faults", default=os.path.join(here, "faults_nwh.json"))
    ap.add_argument("--output", default="tomography_model.xyz")
    ap.add_argument("--figures", default=None, help="directory for QC figures (needs matplotlib)")
    args = ap.parse_args()

    with open(args.faults) as fh:
        spec = json.load(fh)
    settings = spec["settings"]
    faults = [Fault(s, settings["utm_zone"]) for s in spec["faults"]]

    header, data = read_tomography(args.input)
    nums = [l for l in header if l.strip() and not l.lstrip().startswith("#")]
    spacing = [float(v) for v in nums[1].split()]
    dims = [int(v) for v in nums[2].split()]
    print("read %d nodes, grid %s, spacing %s m" % (len(data), dims, spacing))

    reduction, owner = compute_reduction(data, faults, settings, spacing)
    new = apply(data, reduction, settings)
    write_tomography(args.output, header, new)

    print("wrote %s" % args.output)
    print("nodes with Vs reduced: %d (%.2f %%), max reduction %.2f %%"
          % ((reduction > 0).sum(), 100.0 * (reduction > 0).mean(), 100.0 * reduction.max()))
    top = data[:, 2] == data[:, 2].max()
    for k, f in enumerate(faults):
        m = owner == k
        print("  %-5s nodes %7d  (surface nodes at full reduction: %d)"
              % (f.short, m.sum(), (m & top & (reduction >= reduction[top].max() - 1e-9)).sum()))
    print("Vs range: %.2f-%.2f -> %.2f-%.2f m/s"
          % (data[:, 4].min(), data[:, 4].max(), new[:, 4].min(), new[:, 4].max()))
    summarize_geometry(faults, data[:dims[0] * dims[1], :2], settings["end_taper_km"] * 1e3)

    if args.figures:
        make_figures(args.figures, data, new, reduction, faults, dims, settings)
        print("figures in %s" % args.figures)


if __name__ == "__main__":
    main()
