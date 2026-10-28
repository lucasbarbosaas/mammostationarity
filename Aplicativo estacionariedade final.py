"""
MammoStationarity — Análise de Estacionaridade em Mamografia
=============================================================
Instalação dos pacotes usados:
    pip install customtkinter pydicom numpy matplotlib pillow
"""

import threading
import csv
from pathlib import Path
import matplotlib.patches as mpatches

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

import customtkinter as ctk
from tkinter import filedialog, messagebox, Canvas
from PIL import Image, ImageTk
import pydicom

# ==============================================================================
# TEMA
# ==============================================================================

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

C = {
    "bg":    "#0d1117",
    "panel": "#161b22",
    "card":  "#1c2128",
    "border":"#30363d",
    "accent":"#2d7dd2",
    "acc2":  "#1a5fa8",
    "ok":    "#3fb950",
    "text":  "#e6edf3",
    "dim":   "#8b949e",
}

# ==============================================================================
# MOTOR DE ANÁLISE
# ==============================================================================

def load_dicom(filepath):
    ds    = pydicom.dcmread(filepath)
    image = ds.pixel_array.astype(np.float32)
    if hasattr(ds, "RescaleSlope") and hasattr(ds, "RescaleIntercept"):
        image = image * float(ds.RescaleSlope) + float(ds.RescaleIntercept)
    psp = 0.1
    if hasattr(ds, "PixelSpacing"):
        psp = float(ds.PixelSpacing[0])
    elif hasattr(ds, "ImagerPixelSpacing"):
        psp = float(ds.ImagerPixelSpacing[0])
    return image, psp


def get_reference_center(shape, roi_px, psp, ref_dist_mm, chest_wall):
    H, W   = shape
    half   = roi_px // 2
    ref_px = int(round(ref_dist_mm / psp))
    if chest_wall == "left":
        col_ref, row_ref = ref_px, H // 2
    elif chest_wall == "right":
        col_ref, row_ref = W - ref_px, H // 2
    elif chest_wall == "top":
        row_ref, col_ref = ref_px, W // 2
    else:
        row_ref, col_ref = H - ref_px, W // 2
    return (int(np.clip(row_ref, half, H - half - 1)),
            int(np.clip(col_ref, half, W - half - 1)))
def get_reference_center_generic(shape, roi_px):
    H, W = shape
    half = roi_px // 2
    row_ref = int(np.clip(H // 2, half, H - half - 1)) #np.clip para garantir que, em imagens muito pequenas, a ROI não esteja saindo da borda.
    col_ref = int(np.clip(W // 2, half, W - half - 1))
    return row_ref, col_ref

def run_analysis_2d(image, psp, roi_mm, step_px, ref_dist_mm,
                     chest_wall, outlier_pct, progress_cb=None):
    """Varredura 2D completa da imagem."""
    H, W   = image.shape
    roi_px = max(1, int(round(roi_mm / psp)))
    half   = roi_px // 2
    ref_row, ref_col = get_reference_center(
        image.shape, roi_px, psp, ref_dist_mm, chest_wall)

    rows  = np.arange(half, H - half, step_px)
    cols  = np.arange(half, W - half, step_px)
    nr, nc = len(rows), len(cols)

    grid_raw  = np.zeros((nr, nc), dtype=np.float32)
    grid_mean = np.zeros((nr, nc), dtype=np.float32)

    count = 0
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            roi = image[r - half: r + half + 1,
                        c - half: c + half + 1]
            grid_raw[i, j]  = np.std(roi, ddof=0)
            grid_mean[i, j] = np.mean(roi)
            count += 1
            if progress_cb and count % 500 == 0:
                progress_cb(count / (nr * nc))
    if progress_cb:
        progress_cb(1.0)

    # Filtro de outliers
    std_lim  = np.percentile(grid_raw,  outlier_pct)
    mean_lim = np.percentile(grid_mean, outlier_pct)
    mask     = (grid_raw > std_lim) | (grid_mean > mean_lim)
    grid_clean = grid_raw.copy()
    grid_clean[mask] = np.nan

    # Referência
    ri      = int(np.argmin(np.abs(rows - ref_row)))
    ci      = int(np.argmin(np.abs(cols - ref_col)))
    ref_std = float(grid_clean[ri, ci])
    if np.isnan(ref_std) or ref_std == 0:
        ref_std = float(np.nanmedian(grid_clean))

    grid_norm = grid_clean / ref_std
    col_mm    = cols * psp
    row_mm    = (rows - H / 2.0) * psp

    valid = grid_norm[~np.isnan(grid_norm)].ravel()
    stats = {
        "roi_px": roi_px, "roi_mm": roi_mm, "step_px": step_px,
        "n_rois": nr * nc,
        "n_valid": int(np.sum(~np.isnan(grid_norm))),
        "n_removed": int(np.sum(mask)),
        "ref_std": ref_std,
        "min": float(valid.min()), "max": float(valid.max()),
        "mean": float(valid.mean()), "median": float(np.median(valid)),
        "std": float(valid.std()),
        "ratio": float(valid.max() / valid.min()) if valid.min() > 0 else 0,
        "pixel_spacing": psp,
    }
    return grid_norm, grid_raw, row_mm, col_mm, ref_std, stats, rows, cols

def run_analysis_2d_generic(image, psp, roi_mm, step_px, outlier_pct, progress_cb=None):
    """Varredura 2D completa da imagem."""
    H, W   = image.shape
    roi_px = max(1, int(round(roi_mm / psp)))
    half   = roi_px // 2
    ref_row, ref_col = get_reference_center_generic(image.shape, roi_px)

    rows  = np.arange(half, H - half, step_px)
    cols  = np.arange(half, W - half, step_px)
    nr, nc = len(rows), len(cols)

    grid_raw  = np.zeros((nr, nc), dtype=np.float32)
    grid_mean = np.zeros((nr, nc), dtype=np.float32)

    count = 0
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            roi = image[r - half: r + half + 1,
                        c - half: c + half + 1]
            grid_raw[i, j]  = np.std(roi, ddof=0)
            grid_mean[i, j] = np.mean(roi)
            count += 1
            if progress_cb and count % 500 == 0:
                progress_cb(count / (nr * nc))
    if progress_cb:
        progress_cb(1.0)

    # Filtro de outliers
    std_lim  = np.percentile(grid_raw,  outlier_pct)
    mean_lim = np.percentile(grid_mean, outlier_pct)
    mask     = (grid_raw > std_lim) | (grid_mean > mean_lim)
    grid_clean = grid_raw.copy()
    grid_clean[mask] = np.nan

    # Referência
    ri      = int(np.argmin(np.abs(rows - ref_row)))
    ci      = int(np.argmin(np.abs(cols - ref_col)))
    ref_std = float(grid_clean[ri, ci])
    if np.isnan(ref_std) or ref_std == 0:
        ref_std = float(np.nanmedian(grid_clean))

    grid_norm = grid_clean / ref_std
    col_mm    = cols * psp
    row_mm    = (rows - H / 2.0) * psp

    valid = grid_norm[~np.isnan(grid_norm)].ravel()
    stats = {
        "roi_px": roi_px, "roi_mm": roi_mm, "step_px": step_px,
        "n_rois": nr * nc,
        "n_valid": int(np.sum(~np.isnan(grid_norm))),
        "n_removed": int(np.sum(mask)),
        "ref_std": ref_std,
        "min": float(valid.min()), "max": float(valid.max()),
        "mean": float(valid.mean()), "median": float(np.median(valid)),
        "std": float(valid.std()),
        "ratio": float(valid.max() / valid.min()) if valid.min() > 0 else 0,
        "pixel_spacing": psp,
    }
    return grid_norm, grid_raw, row_mm, col_mm, ref_std, stats, rows, cols



def run_line_scan(image, psp, roi_mm, ref_std, chest_wall):
    """
    Linha horizontal na metade exata do eixo Y da imagem,
    varrendo da chest-wall até a borda oposta, passo fixo de 10px. USADA PARA COMPARAÇÃO.
    """
    H, W   = image.shape
    roi_px = max(1, int(round(roi_mm / psp)))
    half   = roi_px // 2
    step   = 10  # px fixo

    # Linha na metade exata do eixo Y
    row_center = H // 2

    # Origem dependendo da chest-wall
    if chest_wall == "left":
        col_start = half
        col_end   = W - half
    elif chest_wall == "right":
        col_start = half
        col_end   = W - half
    elif chest_wall == "top":
        col_start = half
        col_end   = W - half
    else:
        col_start = half
        col_end   = W - half

    scan_cols = np.arange(col_start, col_end, step)

    means = []
    stds  = []
    for c in scan_cols:
        roi = image[row_center - half: row_center + half + 1,
                    c - half:          c + half + 1]
        means.append(float(np.mean(roi)))
        stds.append(float(np.std(roi, ddof=0)))

    means      = np.array(means)
    stds       = np.array(stds)
    stds_norm  = stds / ref_std if ref_std and ref_std > 0 else stds
    col_mm_arr = scan_cols * psp
    row_mm_val = (row_center - H / 2.0) * psp   # deve ser ~0

    return scan_cols, col_mm_arr, means, stds, stds_norm, row_mm_val, row_center

def run_line_scan_generic(image, psp, roi_mm, ref_std):
    """
    Linha horizontal na metade exata do eixo Y da imagem,
    varrendo da chest-wall até a borda oposta, passo fixo de 10px. USADA PARA COMPARAÇÃO.
    """
    H, W   = image.shape
    roi_px = max(1, int(round(roi_mm / psp)))
    half   = roi_px // 2
    step   = 10  # px fixo

    # Linha na metade exata do eixo Y
    row_center = H // 2
    col_start = half
    col_end = W - half
    scan_cols = np.arange(col_start, col_end, step)

    means = []
    stds  = []
    for c in scan_cols:
        roi = image[row_center - half: row_center + half + 1,
                    c - half:          c + half + 1]
        means.append(float(np.mean(roi)))
        stds.append(float(np.std(roi, ddof=0)))

    means      = np.array(means)
    stds       = np.array(stds)
    stds_norm  = stds / ref_std if ref_std and ref_std > 0 else stds
    col_mm_arr = scan_cols * psp
    row_mm_val = (row_center - H / 2.0) * psp   # deve ser ~0

    return scan_cols, col_mm_arr, means, stds, stds_norm, row_mm_val, row_center



# ==============================================================================
# FIGURAS
# ==============================================================================

def make_heatmap_figure(image, grid_norm, row_mm, col_mm, stats, filepath):
    valid  = grid_norm[~np.isnan(grid_norm)].ravel()
    vmin   = np.percentile(valid, 1)
    vmax   = np.percentile(valid, 99)
    masked = np.ma.masked_invalid(grid_norm)
    levels = np.linspace(vmin, vmax, 14)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=100)
    fig.patch.set_facecolor("white")

    ax1 = axes[0]
    ax1.set_facecolor("black")
    ax1.imshow(image, cmap="gray", aspect="equal",
               vmin=np.percentile(image, 1), vmax=np.percentile(image, 99))
    ax1.set_title("Imagem Original", color="white", fontsize=11, pad=8)
    ax1.axis("off")

    ax2 = axes[1]
    ax2.set_facecolor("#7f0000")
    Cm, Rm = np.meshgrid(col_mm, row_mm)
    pm = ax2.pcolormesh(Cm, Rm, masked, cmap="jet",
                         vmin=vmin, vmax=vmax, shading="auto")
    ax2.contour(Cm, Rm, masked, levels=levels,
                colors="k", linewidths=0.3, alpha=0.4)
    # ROI 5×5 mm a 6 cm da parede torácica
    roi_small_w, roi_small_h = 5.0, 5.0
    roi_small_x = 60.0 - roi_small_w / 2
    roi_small_y = -roi_small_h / 2
    ax2.add_patch(mpatches.Rectangle(
        (roi_small_x, roi_small_y), roi_small_w, roi_small_h,
        linewidth=1.5, edgecolor="red", facecolor="none",
        alpha=0.5))

    # ROI 100×100 mm centralizada
    roi_large_w, roi_large_h = 100.0, 100.0
    roi_large_x = col_mm.mean() - roi_large_w / 2
    roi_large_y = -roi_large_h / 2
    ax2.add_patch(mpatches.Rectangle(
        (roi_large_x, roi_large_y), roi_large_w, roi_large_h,
        linewidth=3, edgecolor="black", facecolor="none",
        alpha=0.5))

    # Legenda das ROIs
    leg_roi = ax2.legend(loc="lower right", fontsize=8, framealpha=0.6)
    leg_roi.get_frame().set_facecolor("#1c2128")
    for t in leg_roi.get_texts():
        t.set_color("white")
    #linha abaixo delimita uma área em branco (problema falado na reuniao) onde está a referencia utilizada 6mm parede toraxica
    #if vmin < 1.0 < vmax:
     #   cs = ax2.contour(Cm, Rm, masked, levels=[1.0],
      #                    colors="white", linewidths=1.5, linestyles="--")
       # ax2.clabel(cs, fmt="1.0", fontsize=7, colors="white")

    cbar = fig.colorbar(pm, ax=ax2, fraction=0.046, pad=0.04)
    cbar.set_label("σ / σ_ref", color="black", fontsize=10)
    cbar.ax.yaxis.set_tick_params(color="black")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="black")
    ax2.set_xlabel("Horizontal Position (mm)", color="black", fontsize=9)
    ax2.set_ylabel("Vertical Position (mm)", color="black", fontsize=9)
    ax2.tick_params(colors="black", labelsize=8)
    for sp in ax2.spines.values():
        sp.set_edgecolor("#444")
    ax2.set_aspect("equal")
    ax2.set_title("Stationarity map (σ normalized)",
                  color="black", fontsize=11, pad=8)
    fig.suptitle(
        f"{Path(filepath).name}  |  ROI: {stats['roi_mm']} mm  |  "
        f"Passo: {stats['step_px']} px  |  "
        f"min={stats['min']:.2f}  max={stats['max']:.2f}  razão={stats['ratio']:.2f}",
        color="black", fontsize=9, y=1.01)
    plt.tight_layout()

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    pil = Image.frombuffer("RGBA", canvas.get_width_height(),
                            canvas.buffer_rgba(), "raw", "RGBA", 0, 1)
    return fig, pil


def make_line_figure(col_mm, means, stds, stds_norm, row_mm_val, stats, filepath,
                     xlabel="Posição horizontal (mm)"):
    """Gráfico da varredura horizontal: média e desvio padrão por posição."""
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), dpi=100,
                              gridspec_kw={"hspace": 0.45})
    fig.patch.set_facecolor("#0d1117")

    for ax in axes:
        ax.set_facecolor("white")
        ax.tick_params(colors="black", labelsize=9)
        ax.set_xlabel(xlabel, color="black", fontsize=10)
        for sp in ax.spines.values():
            sp.set_edgecolor("#444")
        ax.grid(color="#30363d", linewidth=0.5, linestyle="--")

    # Painel superior — Média do pixel
    axes[0].plot(col_mm, means, color="#4fc3f7",
                 linewidth=1.4, marker="o", markersize=2.5)
    axes[0].fill_between(col_mm, means, alpha=0.10, color="#4fc3f7")
    axes[0].axhline(np.mean(means), color="white", linestyle="--",
                    linewidth=0.9, alpha=0.6,
                    label=f"Média global: {np.mean(means):.1f}")
    axes[0].set_ylabel("Valor médio do pixel (ADU)", color="white", fontsize=10)
    axes[0].set_title(
        f"Média por ROI  —  linha y = {row_mm_val:.1f} mm (centro da imagem)  "
        f"|  {len(col_mm)} ROIs  |  passo 10 px",
        color="white", fontsize=10)
    leg0 = axes[0].legend(fontsize=9)
    leg0.get_frame().set_facecolor("#1c2128")
    for t in leg0.get_texts():
        t.set_color("white")

    # Painel inferior — σ normalizado
    axes[1].plot(col_mm, stds_norm, color="#ef9a9a",
                 linewidth=1.4, marker="o", markersize=2.5,
                 label="σ / σ_ref")
    axes[1].fill_between(col_mm, stds_norm, alpha=0.10, color="#ef9a9a")
    axes[1].axhline(1.0, color="white", linestyle="--",
                    linewidth=0.9, alpha=0.6, label="Referência (1.0)")
    axes[1].axhline(np.mean(stds_norm), color="#ffd54f", linestyle=":",
                    linewidth=1.0,
                    label=f"Média: {np.mean(stds_norm):.3f}")
    axes[1].set_ylabel("σ normalizado  (σ / σ_ref)", color="white", fontsize=10)
    axes[1].set_title("Desvio Padrão normalizado por ROI", color="white", fontsize=10)
    leg1 = axes[1].legend(fontsize=9)
    leg1.get_frame().set_facecolor("#1c2128")
    for t in leg1.get_texts():
        t.set_color("white")

    fig.suptitle(
        f"{Path(filepath).name}  |  ROI: {stats['roi_mm']} mm  |  "
        f"σ_ref = {stats['ref_std']:.4f}",
        color="#8b949e", fontsize=9)
    plt.tight_layout()

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    pil = Image.frombuffer("RGBA", canvas.get_width_height(),
                            canvas.buffer_rgba(), "raw", "RGBA", 0, 1)
    return fig, pil


# ==============================================================================
# WIDGET ZOOM / PAN
# ==============================================================================

class ZoomCanvas(ctk.CTkFrame):
    def __init__(self, master, **kwargs):
        super().__init__(master, fg_color=C["card"], corner_radius=8, **kwargs)
        self._cv = Canvas(self, bg=C["card"], highlightthickness=0, cursor="fleur")
        self._cv.pack(fill="both", expand=True)

        bar = ctk.CTkFrame(self, fg_color=C["panel"], corner_radius=0, height=30)
        bar.pack(fill="x", side="bottom")
        bkw = dict(width=30, height=24, fg_color=C["card"],
                   hover_color=C["border"], border_width=1,
                   border_color=C["border"], text_color=C["text"],
                   font=ctk.CTkFont(size=13))
        ctk.CTkButton(bar, text="＋", command=self._zi, **bkw).pack(side="left", padx=4, pady=3)
        ctk.CTkButton(bar, text="－", command=self._zo, **bkw).pack(side="left", pady=3)
        ctk.CTkButton(bar, text="⊡",  command=self._fit,**bkw).pack(side="left", padx=4, pady=3)
        self._lbl = ctk.CTkLabel(bar, text="100%", font=ctk.CTkFont(size=10),
                                  text_color=C["dim"])
        self._lbl.pack(side="left", padx=6)

        self._img = None; self._tk = None
        self._sc  = 1.0;  self._ox = 0; self._oy = 0; self._drag = None

        self._cv.bind("<ButtonPress-1>", lambda e: setattr(self, "_drag", (e.x, e.y)))
        self._cv.bind("<B1-Motion>",     self._pan)
        self._cv.bind("<MouseWheel>",    self._scroll)
        self._cv.bind("<Button-4>",      self._scroll)
        self._cv.bind("<Button-5>",      self._scroll)
        self.bind("<Configure>",         lambda e: self._fit())

    def set_image(self, pil_img):
        self._img = pil_img
        self.after(60, self._fit)

    def _fit(self):
        if self._img is None: return
        cw, ch = self._cv.winfo_width(), self._cv.winfo_height()
        if cw < 10 or ch < 10:
            self.after(100, self._fit); return
        iw, ih = self._img.size
        self._sc = min(cw / iw, ch / ih)
        self._ox = (cw - iw * self._sc) / 2
        self._oy = (ch - ih * self._sc) / 2
        self._render()

    def _render(self):
        if self._img is None: return
        iw, ih = self._img.size
        nw, nh = max(1, int(iw * self._sc)), max(1, int(ih * self._sc))
        self._tk = ImageTk.PhotoImage(self._img.resize((nw, nh), Image.LANCZOS))
        self._cv.delete("all")
        self._cv.create_image(self._ox, self._oy, anchor="nw", image=self._tk)
        self._lbl.configure(text=f"{self._sc*100:.0f}%")

    def _zoom(self, f, cx=None, cy=None):
        if self._img is None: return
        cx = cx or self._cv.winfo_width() / 2
        cy = cy or self._cv.winfo_height() / 2
        ns = max(0.05, min(self._sc * f, 20.0))
        r  = ns / self._sc
        self._ox = cx - r * (cx - self._ox)
        self._oy = cy - r * (cy - self._oy)
        self._sc = ns
        self._render()

    def _zi(self):     self._zoom(1.25)
    def _zo(self):     self._zoom(0.8)
    def _scroll(self, e):
        self._zoom(1.15 if (e.num == 4 or e.delta > 0) else 1/1.15, e.x, e.y)
    def _pan(self, e):
        if not self._drag: return
        self._ox += e.x - self._drag[0]
        self._oy += e.y - self._drag[1]
        self._drag = (e.x, e.y)
        self._render()


# ==============================================================================
# PAINEL ESQUERDO COM SCROLL
# ==============================================================================

class ScrollableLeft(ctk.CTkScrollableFrame):
    """Painel de controles com scroll vertical automático."""
    def __init__(self, master, app, **kwargs):
        super().__init__(master, fg_color=C["panel"],
                          scrollbar_button_color=C["border"],
                          scrollbar_button_hover_color=C["accent"],
                          **kwargs)
        self._app = app
        self._build()

    def _section(self, title):
        ctk.CTkLabel(self, text=title.upper(),
                     font=ctk.CTkFont(size=10, weight="bold"),
                     text_color=C["dim"]).pack(anchor="w", padx=14, pady=(12, 1))
        ctk.CTkFrame(self, height=1, fg_color=C["border"]).pack(
            fill="x", padx=14, pady=(0, 6))

    def _row_entry(self, label, default, width=78, master=None):
        parent = master if master is not None else self
        f = ctk.CTkFrame(parent, fg_color="transparent")
        f.pack(fill="x", padx=14, pady=2)
        ctk.CTkLabel(f, text=label, font=ctk.CTkFont(size=12),
                     text_color=C["text"], width=148, anchor="w").pack(side="left")
        var = ctk.StringVar(value=str(default))
        ctk.CTkEntry(f, textvariable=var, width=width,
                     fg_color=C["card"], border_color=C["border"],
                     text_color=C["text"]).pack(side="right")
        return var

    def _build(self):
        app = self._app

        # Arquivo
        self._section("Arquivo DICOM")
        app.lbl_file = ctk.CTkLabel(
            self, text="Nenhum arquivo selecionado",
            font=ctk.CTkFont(size=11), text_color=C["dim"],
            wraplength=220, justify="left")
        app.lbl_file.pack(anchor="w", padx=14, pady=(0, 5))
        ctk.CTkButton(self, text="📂  Selecionar DICOM",
                      command=app._select_file,
                      fg_color=C["accent"], hover_color=C["acc2"],
                      font=ctk.CTkFont(size=12, weight="bold"),
                      height=36).pack(fill="x", padx=14, pady=(0, 4))

        #Modo de Análise
        self._section("Modo de Análise")

        toggle_frame = ctk.CTkFrame(self, fg_color="transparent")
        toggle_frame.pack(fill="x", padx=14, pady=(0, 6))

        app.btn_mode_mammo = ctk.CTkButton(
            toggle_frame, text="Mamografia",
            command=lambda: app._set_mode("mammo"),
            fg_color=C["accent"], hover_color=C["acc2"],
            font=ctk.CTkFont(size=12), height=30, width=116)
        app.btn_mode_mammo.pack(side="left")

        app.btn_mode_generic = ctk.CTkButton(
            toggle_frame, text="Genérico",
            command=lambda: app._set_mode("generic"),
            fg_color=C["card"], hover_color=C["border"],
            border_width=1, border_color=C["border"],
            font=ctk.CTkFont(size=12), height=30, width=116)
        app.btn_mode_generic.pack(side="right")

        # Parâmetros
        self._section("Parâmetros")
        app.var_step = self._row_entry("Passo 2D (px)", 10)
        app.var_out = self._row_entry("Filtro outlier (%)", 97)
        app.frame_mammo_params = ctk.CTkFrame(self, fg_color="transparent")
        app.frame_mammo_params.pack(fill="x")
        app.var_roi = self._row_entry("ROI (mm)", 5.0, master=app.frame_mammo_params)
        app.var_ref = self._row_entry("Dist. referência (mm)", 60, master=app.frame_mammo_params)
        app.frame_generic_params = ctk.CTkFrame(self, fg_color="transparent")

        app.var_roi_generic = self._row_entry("ROI (mm)", 10.0, master=app.frame_generic_params)

        f_cw = ctk.CTkFrame(app.frame_mammo_params, fg_color="transparent")
        f_cw.pack(fill="x", padx=14, pady=2)
        ctk.CTkLabel(f_cw, text="Chest-wall",
                     font=ctk.CTkFont(size=12), text_color=C["text"],
                     width=148, anchor="w").pack(side="left")
        app.var_cw = ctk.StringVar(value="left")
        ctk.CTkOptionMenu(f_cw, values=["left", "right", "top", "bottom"],
                          variable=app.var_cw, width=78,
                          fg_color=C["card"], button_color=C["accent"],
                          button_hover_color=C["acc2"],
                          text_color=C["text"]).pack(side="right")

        # Analisar
        self._section("Análise")
        app.btn_run = ctk.CTkButton(
            self, text="▶  Analisar", command=app._start,
            fg_color=C["ok"], hover_color="#2ea043",
            font=ctk.CTkFont(size=13, weight="bold"), height=40)
        app.btn_run.pack(fill="x", padx=14, pady=(0, 6))

        app.progress = ctk.CTkProgressBar(self, height=7,
                                           fg_color=C["card"],
                                           progress_color=C["accent"])
        app.progress.pack(fill="x", padx=14, pady=(0, 2))
        app.progress.set(0)

        app.lbl_prog = ctk.CTkLabel(self, text="",
                                     font=ctk.CTkFont(size=10),
                                     text_color=C["dim"])
        app.lbl_prog.pack(anchor="w", padx=14)

        # Estatísticas
        self._section("Resultados")
        sf = ctk.CTkFrame(self, fg_color=C["card"], corner_radius=7)
        sf.pack(fill="x", padx=14, pady=(0, 6))
        app.lbl_stats = ctk.CTkLabel(
            sf, text="—", font=ctk.CTkFont("Courier New", 11),
            text_color=C["dim"], justify="left")
        app.lbl_stats.pack(padx=10, pady=8, anchor="w")

        # Exportar
        self._section("Exportar")

        def export_btn(text, cmd):
            b = ctk.CTkButton(
                self, text=text, command=cmd,
                fg_color=C["card"], hover_color=C["border"],
                border_width=1, border_color=C["accent"],
                text_color=C["accent"],
                font=ctk.CTkFont(size=12), height=34, state="disabled")
            b.pack(fill="x", padx=14, pady=(0, 5))
            return b

        app.btn_png      = export_btn("💾  Salvar PNG (mapa 2D)",   app._save_png)
        app.btn_png_line = export_btn("💾  Salvar PNG (linha)",     app._save_png_line)
        app.btn_csv      = export_btn("📊  Salvar CSV (mapa 2D)",   app._save_csv)
        app.btn_csv_line = export_btn("📊  Salvar CSV (linha)",     app._save_csv_line)

        # Espaço extra no final para o scroll não cortar o último botão
        ctk.CTkLabel(self, text="", height=20).pack()


# ==============================================================================
# APLICATIVO PRINCIPAL
# ==============================================================================

class MammoApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("CalculadoraEstacionaridade")
        self.geometry("1280x720")
        self.minsize(1000, 660)
        self.configure(fg_color=C["bg"])

        self.filepath   = None
        self.image      = None
        self.psp        = None
        self.grid_norm  = None
        self.grid_raw   = None
        self.row_mm     = None
        self.col_mm     = None
        self.rows       = None
        self.cols       = None
        self.stats      = None
        self.fig_heat   = None
        self.fig_line   = None
        self.line_data  = None  # dict com arrays da linha
        self._running   = False
        self.mode = "mammo"  # "mammo" ou "generic"
        self.frame_mammo_params = None  # frame com params exclusivos da mamografia
        self.frame_generic_params = None  # frame com params exclusivos do modo genérico
        self.var_roi_generic = None  # ROI separado para o modo genérico

        # Widgets criados pelo ScrollableLeft
        self.lbl_file   = None
        self.var_roi    = None
        self.var_step   = None
        self.var_ref    = None
        self.var_out    = None
        self.var_cw     = None
        self.btn_run    = None
        self.progress   = None
        self.lbl_prog   = None
        self.lbl_stats  = None
        self.btn_png    = None
        self.btn_png_line = None
        self.btn_csv    = None
        self.btn_csv_line = None

        self._build_ui()

    def _set_mode(self, mode):
        self.mode = mode
        if mode == "mammo":
            self.frame_generic_params.pack_forget()
            self.frame_mammo_params.pack(fill="x")
            self.btn_mode_mammo.configure(fg_color=C["accent"])
            self.btn_mode_generic.configure(fg_color=C["card"])
        else:
            self.frame_mammo_params.pack_forget()
            self.frame_generic_params.pack(fill="x")
            self.btn_mode_mammo.configure(fg_color=C["card"])
            self.btn_mode_generic.configure(fg_color=C["accent"])

    def _build_ui(self):
        # Cabeçalho
        hdr = ctk.CTkFrame(self, fg_color=C["panel"], corner_radius=0, height=52)
        hdr.pack(fill="x")
        hdr.pack_propagate(False)
        ctk.CTkLabel(hdr, text="⬡  CalculadoraEstacionaridade",
                     font=ctk.CTkFont("Courier New", 17, "bold"),
                     text_color=C["accent"]).pack(side="left", padx=18, pady=8)
        ctk.CTkLabel(hdr,
                     text="Análise de Estacionaridade Espacial de Ruído — DICOM",
                     font=ctk.CTkFont(size=11),
                     text_color=C["dim"]).pack(side="left", padx=2)

        body = ctk.CTkFrame(self, fg_color="transparent")
        body.pack(fill="both", expand=True, padx=10, pady=8)
        body.columnconfigure(1, weight=1)
        body.rowconfigure(0, weight=1)

        # Painel esquerdo com scroll
        ScrollableLeft(body, app=self, width=272).grid(
            row=0, column=0, sticky="nsew", padx=(0, 8))

        # Painel direito com abas
        right = ctk.CTkFrame(body, fg_color=C["panel"], corner_radius=10)
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)
        self._build_tabs(right)

    def _build_tabs(self, parent):
        tab_bar = ctk.CTkFrame(parent, fg_color=C["card"],
                                corner_radius=0, height=38)
        tab_bar.grid(row=0, column=0, sticky="ew", padx=6, pady=(6, 0))
        tab_bar.pack_propagate(False)

        self.btn_tab_heat = ctk.CTkButton(
            tab_bar, text="Mapa de Calor",
            command=lambda: self._switch_tab("heat"),
            fg_color=C["accent"], hover_color=C["acc2"],
            font=ctk.CTkFont(size=12), height=30, width=160)
        self.btn_tab_heat.pack(side="left", padx=6, pady=4)

        self.btn_tab_line = ctk.CTkButton(
            tab_bar, text="Linha Horizontal",
            command=lambda: self._switch_tab("line"),
            fg_color=C["card"], hover_color=C["border"],
            font=ctk.CTkFont(size=12), height=30, width=160)
        self.btn_tab_line.pack(side="left", pady=4)

        tc = ctk.CTkFrame(parent, fg_color="transparent")
        tc.grid(row=1, column=0, sticky="nsew", padx=6, pady=6)
        tc.rowconfigure(0, weight=1)
        tc.columnconfigure(0, weight=1)
        self.tab_container = tc

        # Aba mapa
        self.tab_heat = ctk.CTkFrame(tc, fg_color="transparent")
        self.tab_heat.grid(row=0, column=0, sticky="nsew")
        self.tab_heat.rowconfigure(0, weight=1)
        self.tab_heat.columnconfigure(0, weight=1)
        self.ph_heat = ctk.CTkLabel(
            self.tab_heat,
            text="Selecione um arquivo DICOM\ne clique em  ▶ Analisar",
            font=ctk.CTkFont(size=14), text_color=C["dim"])
        self.ph_heat.grid(row=0, column=0)
        self.zc_heat = ZoomCanvas(self.tab_heat)

        # Aba linha
        self.tab_line = ctk.CTkFrame(tc, fg_color="transparent")
        self.tab_line.rowconfigure(0, weight=1)
        self.tab_line.columnconfigure(0, weight=1)
        self.ph_line = ctk.CTkLabel(
            self.tab_line,
            text="Execute a análise para ver\no gráfico da linha horizontal",
            font=ctk.CTkFont(size=14), text_color=C["dim"])
        self.ph_line.grid(row=0, column=0)
        self.zc_line = ZoomCanvas(self.tab_line)

        self._cur_tab = "heat"

    def _switch_tab(self, tab):
        self._cur_tab = tab
        if tab == "heat":
            self.tab_line.grid_forget()
            self.tab_heat.grid(row=0, column=0, sticky="nsew")
            self.btn_tab_heat.configure(fg_color=C["accent"])
            self.btn_tab_line.configure(fg_color=C["card"])
        else:
            self.tab_heat.grid_forget()
            self.tab_line.grid(row=0, column=0, sticky="nsew")
            self.btn_tab_line.configure(fg_color=C["accent"])
            self.btn_tab_heat.configure(fg_color=C["card"])

    # ── Lógica ────────────────────────────────────────────────────────

    def _select_file(self):
        path = filedialog.askopenfilename(
            title="Selecionar arquivo DICOM",
            filetypes=[("Todos os arquivos", "*.*"),
                       ("DICOM files", "*.dcm *.DCM")])
        if not path:
            return
        self.filepath = path
        self.lbl_file.configure(
            text=Path(path).name or Path(path).parts[-1],
            text_color=C["text"])
        self.progress.set(0)
        self.lbl_prog.configure(text="")

    def _start(self):
        if not self.filepath:
            messagebox.showwarning("Atenção", "Selecione um arquivo DICOM primeiro.")
            return
        if self._running:
            return

        try:
            step_px = int(self.var_step.get())  # compartilhado
            outlier = float(self.var_out.get())  # compartilhado

            if self.mode == "mammo":
                roi_mm = float(self.var_roi.get())
                ref_dist = float(self.var_ref.get())
                chest_wall = self.var_cw.get()
            else:
                roi_mm = float(self.var_roi_generic.get())
                # ref_dist e chest_wall não existem aqui

        except ValueError:
            messagebox.showerror("Erro", "Verifique os parâmetros.")
            return

        self._running = True
        self.btn_run.configure(state="disabled", text="Analisando...")
        for b in [self.btn_png, self.btn_png_line,
                  self.btn_csv, self.btn_csv_line]:
            b.configure(state="disabled")
        self.progress.set(0)
        self.lbl_prog.configure(text="Carregando DICOM...")

        def worker():
            try:
                # PRIMEIRO: carrega a imagem — isso é comum aos dois modos
                image, psp = load_dicom(self.filepath)
                self.image = image
                self.psp = psp

                def on_prog(p):
                    self.after(0, lambda v=p: self.progress.set(v))
                    self.after(0, lambda v=p: self.lbl_prog.configure(
                        text=f"{v * 100:.0f}%"))

                # DEPOIS: bifurca dependendo do modo
                if self.mode == "mammo":
                    (self.grid_norm, self.grid_raw,
                     self.row_mm, self.col_mm,
                     ref_std, self.stats,
                     self.rows, self.cols) = run_analysis_2d(
                        image, psp, roi_mm, step_px,
                        ref_dist, chest_wall, outlier,
                        progress_cb=on_prog)
                    scan_cols, col_mm_line, means, stds, stds_norm, row_mm_val, row_center = \
                        run_line_scan(image, psp, roi_mm, ref_std, chest_wall)
                    xlabel = "Posição horizontal — da chest-wall (mm)"
                else:
                    (self.grid_norm, self.grid_raw,
                     self.row_mm, self.col_mm,
                     ref_std, self.stats,
                     self.rows, self.cols) = run_analysis_2d_generic(
                        image, psp, roi_mm, step_px,
                        outlier,
                        progress_cb=on_prog)
                    scan_cols, col_mm_line, means, stds, stds_norm, row_mm_val, row_center = \
                        run_line_scan_generic(image, psp, roi_mm, ref_std)
                    xlabel = "Posição horizontal (mm)"

                # DEPOIS da bifurcação: figuras e atualização da UI
                self.line_data = {
                    "scan_cols": scan_cols,
                    "col_mm": col_mm_line,
                    "means": means,
                    "stds": stds,
                    "stds_norm": stds_norm,
                    "row_mm_val": row_mm_val,
                    "row_center": row_center,
                }
                fig_h, pil_h = make_heatmap_figure(
                    image, self.grid_norm, self.row_mm, self.col_mm,
                    self.stats, self.filepath)
                self.fig_heat = fig_h
                fig_l, pil_l = make_line_figure(
                    col_mm_line, means, stds, stds_norm,
                    row_mm_val, self.stats, self.filepath,
                    xlabel=xlabel)
                self.fig_line = fig_l
                self.after(0, lambda: self._update(pil_h, pil_l))

            except Exception as e:
                self.after(0, lambda err=e: messagebox.showerror("Erro", str(err)))
            finally:
                self._running = False
                self.after(0, lambda: self.btn_run.configure(
                    state="normal", text="▶  Analisar"))

        threading.Thread(target=worker, daemon=True).start()

    def _update(self, pil_heat, pil_line):
        s = self.stats
        txt = (
            f"ROIs totais  : {s['n_rois']}\n"
            f"ROIs válidas : {s['n_valid']}\n"
            f"Removidas    : {s['n_removed']}\n"
            f"─────────────────\n"
            f"σ_ref        : {s['ref_std']:.4f}\n"
            f"mín          : {s['min']:.4f}\n"
            f"máx          : {s['max']:.4f}\n"
            f"média        : {s['mean']:.4f}\n"
            f"mediana      : {s['median']:.4f}\n"
            f"razão        : {s['ratio']:.3f}\n"
            f"─────────────────\n"
            f"px spacing   : {s['pixel_spacing']:.4f} mm"
        )
        self.lbl_stats.configure(text=txt, text_color=C["text"])
        self.lbl_prog.configure(text="✓  Concluído", text_color=C["ok"])

        self.ph_heat.grid_forget()
        self.zc_heat.grid(row=0, column=0, sticky="nsew")
        self.zc_heat.set_image(pil_heat)

        self.ph_line.grid_forget()
        self.zc_line.grid(row=0, column=0, sticky="nsew")
        self.zc_line.set_image(pil_line)

        for b in [self.btn_png, self.btn_png_line,
                  self.btn_csv, self.btn_csv_line]:
            b.configure(state="normal")

    # ── Exportar ──────────────────────────────────────────────────────

    def _save_png(self):
        if not self.fig_heat: return
        path = filedialog.asksaveasfilename(
            title="Salvar mapa de calor",
            defaultextension=".png",
            initialfile=Path(self.filepath).stem + "_stationarity_map.png",
            filetypes=[("PNG", "*.png"), ("Todos", "*.*")])
        if path:
            self.fig_heat.savefig(path, dpi=200, bbox_inches="tight",
                                   facecolor=self.fig_heat.get_facecolor())
            messagebox.showinfo("Salvo", f"PNG salvo:\n{path}")

    def _save_png_line(self):
        if not self.fig_line: return
        path = filedialog.asksaveasfilename(
            title="Salvar gráfico da linha",
            defaultextension=".png",
            initialfile=Path(self.filepath).stem + "_linha_horizontal.png",
            filetypes=[("PNG", "*.png"), ("Todos", "*.*")])
        if path:
            self.fig_line.savefig(path, dpi=200, bbox_inches="tight",
                                   facecolor=self.fig_line.get_facecolor())
            messagebox.showinfo("Salvo", f"PNG salvo:\n{path}")

    def _save_csv(self):
        """CSV 2D completo."""
        if self.grid_norm is None: return
        path = filedialog.asksaveasfilename(
            title="Salvar CSV do mapa 2D",
            defaultextension=".csv",
            initialfile=Path(self.filepath).stem + "_stationarity_2D.csv",
            filetypes=[("CSV", "*.csv"), ("Todos", "*.*")])
        if not path: return
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["row_px", "col_px", "row_mm", "col_mm",
                         "std_raw", "std_normalizado", "outlier_removido"])
            for i, r in enumerate(self.rows):
                for j, c in enumerate(self.cols):
                    raw  = self.grid_raw[i, j]
                    norm = self.grid_norm[i, j]
                    w.writerow([int(r), int(c),
                                f"{self.row_mm[i]:.3f}",
                                f"{self.col_mm[j]:.3f}",
                                f"{raw:.6f}",
                                f"{norm:.6f}" if not np.isnan(norm) else "",
                                "sim" if np.isnan(norm) else "nao"])
        messagebox.showinfo("Salvo", f"CSV salvo:\n{path}")

    def _save_csv_line(self):
        """CSV da linha horizontal — pronto para comparar com ImageJ."""
        if not self.line_data: return
        ld   = self.line_data
        path = filedialog.asksaveasfilename(
            title="Salvar CSV da linha horizontal",
            defaultextension=".csv",
            initialfile=Path(self.filepath).stem + "_linha_horizontal.csv",
            filetypes=[("CSV", "*.csv"), ("Todos", "*.*")])
        if not path: return
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["col_px", "col_mm",
                         "media_pixel", "std_raw", "std_normalizado"])
            for col_px, mm, mean, std, norm in zip(
                    ld["scan_cols"], ld["col_mm"],
                    ld["means"], ld["stds"], ld["stds_norm"]):
                w.writerow([int(col_px), f"{mm:.3f}",
                             f"{mean:.4f}", f"{std:.6f}", f"{norm:.6f}"])
        messagebox.showinfo("Salvo",
            f"CSV da linha salvo:\n{path}\n\n"
            f"Linha em y = {ld['row_mm_val']:.1f} mm  "
            f"(pixel row = {ld['row_center']})\n"
            f"{len(ld['col_mm'])} ROIs  |  passo 10 px")


# ==============================================================================
# MAIN
# ==============================================================================

if __name__ == "__main__":
    app = MammoApp()
    app.mainloop()
