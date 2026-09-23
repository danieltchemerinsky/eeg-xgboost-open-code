# ==============================================================================
# FIGURE 1 - PATIENT INCLUSION FLOWCHART
# ==============================================================================
# Draws the flowchart, deriving its counts from the data wherever they can be
# derived rather than declaring them.
#
# Of the six counts in the figure, four come from outputminus.RData:
#
#   patients with EEG at 12-24 hours   nrow(combined)
#   not included in the ML model       those with > 5 missing predictors
#   patients included in the ML model  the remainder
#   good / poor outcome                CPC 1-3 vs CPC 4-5 among those included
#
# The two upstream counts - patients recruited at the two sites, and those
# missed in the EEG substudy - are not recoverable from this dataset, which
# begins after that exclusion. They are declared as constants below and are the
# only numbers in the figure that must be checked against the study records.
# The script verifies that they are consistent with the data.
#
# Every element is placed at an explicit millimetre coordinate rather than by an
# automatic layout engine, so box widths, column edges and the split under the
# ML box align exactly and stay aligned if the text changes.
#
# Run standalone (it does not need 01 to have been run):
#   Rscript 09_figure1_flowchart.R
#
# Uses only `grid` and `grDevices`, both part of base R.
# ==============================================================================

library(grid)

# ------------------------------------------------------------------------------
# Counts that come from the study records, not from the data file
# ------------------------------------------------------------------------------
N_RECRUITED    <- 159
N_MISSED       <- 38
MISSED_REASONS <- c("Late inclusion in TH48",
                    "Early transfer of patients",
                    "Logistical reasons")

DATA_FILE   <- "outputminus.RData"
MAX_MISSING <- 5          # same threshold as 01_loocv_analysis.R
OUT_STEM    <- "Figure1_Flowchart"

# ------------------------------------------------------------------------------
# Derive the counts from the data
# ------------------------------------------------------------------------------
if (!file.exists(DATA_FILE))
  stop(DATA_FILE, " not found in ", getwd())

.e <- new.env(); load(DATA_FILE, envir = .e); combined <- .e$combined; rm(.e)

n_eeg    <- nrow(combined)
keep     <- rowSums(is.na(combined)) <= MAX_MISSING
n_excl   <- sum(!keep)
n_model  <- sum(keep)
outcomes <- combined$Outcome[keep]
n_good   <- sum(outcomes %in% c(1, 2, 3))
n_poor   <- sum(outcomes %in% c(4, 5))

cat("Derived from ", DATA_FILE, ":\n", sep = "")
cat(sprintf("  EEG at 12-24 hours     : %d\n", n_eeg))
cat(sprintf("  > %d missing predictors : %d\n", MAX_MISSING, n_excl))
cat(sprintf("  included in ML model   : %d\n", n_model))
cat(sprintf("  good (CPC 1-3)         : %d\n", n_good))
cat(sprintf("  poor (CPC 4-5)         : %d\n", n_poor))
cat("From the study records:\n")
cat(sprintf("  recruited at two sites : %d\n", N_RECRUITED))
cat(sprintf("  missed in EEG substudy : %d\n\n", N_MISSED))

if (N_RECRUITED - N_MISSED != n_eeg)
  warning(sprintf("%d recruited minus %d missed is %d, but the data hold %d.",
                  N_RECRUITED, N_MISSED, N_RECRUITED - N_MISSED, n_eeg))
if (n_good + n_poor != n_model)
  warning("Good and poor outcomes do not sum to the number included.")

# ==============================================================================
# LAYOUT - all coordinates in millimetres on a 195 x 124 mm canvas
# ==============================================================================
CANVAS_W <- 195; CANVAS_H <- 124

CX <- 58;  CW <- 88          # main column centre and width
SX <- 149; SW <- 78          # annotation column centre and width
BH <- 17                     # standard box height
X1H <- 26                    # taller annotation box

Y_A <- 105; Y_B <- 78; Y_C <- 51; Y_O <- 15; Y_CONN <- 30

OUT_W <- 42; OUT_GAP <- 4
GX <- CX - OUT_GAP/2 - OUT_W/2
PX <- CX + OUT_GAP/2 + OUT_W/2

BLUE  <- list(fill = "#DCEBFA", line = "#2E6FB7", text = "#14365C")
TEAL  <- list(fill = "#D6F2EC", line = "#2E8B87", text = "#123F3D")
GREEN <- list(fill = "#C3E6C6", line = "#3E9E44", text = "#1B4620")
RED   <- list(fill = "#F6C9C9", line = "#D7403B", text = "#5E1714")
NOTE  <- list(fill = "#FDF3F3", line = "#D9534F", text = "#8E211D")
ARROW_COL <- "#31445B"; DASH_COL <- "#D9534F"

mm <- function(v) unit(v, "native")

box <- function(cx, cy, w, h, style, dashed = FALSE) {
  grid.roundrect(x = mm(cx), y = mm(cy), width = mm(w), height = mm(h),
                 r = unit(3, "mm"),
                 gp = gpar(fill = style$fill, col = style$line, lwd = 1.5,
                           lty = if (dashed) "dashed" else "solid"))
}

two_line <- function(cx, cy, title, sub, style, ts = 10.5, ss = 9) {
  grid.text(title, x = mm(cx), y = mm(cy + 3.0),
            gp = gpar(fontsize = ts, fontface = "bold", col = style$text,
                      fontfamily = "sans"))
  grid.text(sub, x = mm(cx), y = mm(cy - 3.4),
            gp = gpar(fontsize = ss, col = style$text, fontfamily = "sans"))
}

# Arrowheads are drawn in absolute units so they stay the same size whatever
# the canvas is scaled to.
HEAD <- arrow(angle = 22, length = unit(2.4, "mm"), ends = "last", type = "closed")

arrow_v <- function(x, y0, y1, col = ARROW_COL) {
  grid.lines(x = mm(c(x, x)), y = mm(c(y0, y1)),
             arrow = HEAD, gp = gpar(col = col, fill = col, lwd = 1.5))
}

plain_line <- function(x0, y0, x1, y1, col = ARROW_COL, lty = "solid") {
  grid.lines(x = mm(c(x0, x1)), y = mm(c(y0, y1)),
             gp = gpar(col = col, lwd = 1.5, lty = lty, lineend = "round"))
}

# A dashed shaft with a solid head: drawing the whole thing with lty = "dashed"
# would break up the arrowhead outline too.
dashed_arrow <- function(x0, x1, y) {
  plain_line(x0, y, x1 - 3, y, col = DASH_COL, lty = "22")
  grid.lines(x = mm(c(x1 - 3, x1)), y = mm(c(y, y)),
             arrow = HEAD, gp = gpar(col = DASH_COL, fill = DASH_COL, lwd = 1.5))
}

draw_flowchart <- function() {
  grid.newpage()
  pushViewport(viewport(width = unit(CANVAS_W, "mm"), height = unit(CANVAS_H, "mm"),
                        xscale = c(0, CANVAS_W), yscale = c(0, CANVAS_H)))
  grid.rect(gp = gpar(fill = "white", col = NA))

  # --- connectors first, so the boxes sit on top of the line ends ------------
  arrow_v(CX, Y_A - BH/2, Y_B + BH/2)
  arrow_v(CX, Y_B - BH/2, Y_C + BH/2)
  plain_line(CX, Y_C - BH/2, CX, Y_CONN)
  plain_line(GX, Y_CONN, PX, Y_CONN)
  arrow_v(GX, Y_CONN, Y_O + BH/2)
  arrow_v(PX, Y_CONN, Y_O + BH/2)
  dashed_arrow(CX + CW/2, SX - SW/2, Y_A)
  dashed_arrow(CX + CW/2, SX - SW/2, Y_B)

  # --- main column ------------------------------------------------------------
  box(CX, Y_A, CW, BH, BLUE)
  two_line(CX, Y_A, "Patients included at two sites",
           sprintf("(n=%d)", N_RECRUITED), BLUE)

  box(CX, Y_B, CW, BH, BLUE)
  two_line(CX, Y_B, "Patients with EEG at 12-24 hours",
           sprintf("(n=%d)", n_eeg), BLUE)

  box(CX, Y_C, CW, BH, TEAL)
  two_line(CX, Y_C, "Patients included in ML model",
           sprintf("(n=%d)", n_model), TEAL)

  # --- outcomes ---------------------------------------------------------------
  box(GX, Y_O, OUT_W, BH, GREEN)
  two_line(GX, Y_O, "Good outcome",
           sprintf("CPC 1-3 (n=%d)", n_good), GREEN)

  box(PX, Y_O, OUT_W, BH, RED)
  two_line(PX, Y_O, "Poor outcome",
           sprintf("CPC 4-5 (n=%d)", n_poor), RED)

  # --- annotations ------------------------------------------------------------
  box(SX, Y_A, SW, X1H, NOTE, dashed = TRUE)
  grid.text(sprintf("Missed in EEG substudy (n=%d)", N_MISSED),
            x = mm(SX), y = mm(Y_A + 8.4),
            gp = gpar(fontsize = 9.5, fontface = "bold", col = NOTE$text,
                      fontfamily = "sans"))
  # Centred, matching the box below. Text is ASCII throughout: the default
  # pdf() device silently substitutes characters its font lacks.
  for (i in seq_along(MISSED_REASONS))
    grid.text(MISSED_REASONS[i],
              x = mm(SX), y = mm(Y_A + 2.2 - (i - 1) * 5.0),
              gp = gpar(fontsize = 8.8, col = NOTE$text, fontfamily = "sans"))

  box(SX, Y_B, SW, BH, NOTE, dashed = TRUE)
  grid.text(sprintf("Not included in ML model (n=%d)", n_excl),
            x = mm(SX), y = mm(Y_B + 3.0),
            gp = gpar(fontsize = 9.5, fontface = "bold", col = NOTE$text,
                      fontfamily = "sans"))
  grid.text(sprintf("More than %d missing predictors", MAX_MISSING),
            x = mm(SX), y = mm(Y_B - 3.4),
            gp = gpar(fontsize = 8.8, col = NOTE$text, fontfamily = "sans"))

  popViewport()
}

# ==============================================================================
# RENDER
# ==============================================================================
w_in <- CANVAS_W / 25.4; h_in <- CANVAS_H / 25.4

pdf(paste0(OUT_STEM, ".pdf"), width = w_in, height = h_in)
draw_flowchart(); invisible(dev.off())
cat("Saved: ", OUT_STEM, ".pdf\n", sep = "")

# The device default is used deliberately: asking for type = "cairo" fails on a
# macOS install without XQuartz, and the fallback then loses glyphs.
png(paste0(OUT_STEM, ".png"), width = CANVAS_W, height = CANVAS_H,
    units = "mm", res = 600)
draw_flowchart(); invisible(dev.off())
cat("Saved: ", OUT_STEM, ".png\n", sep = "")

# SVG needs cairo, or the svglite package. Skip quietly if neither is present -
# the PDF is already vector, and is what a journal will ask for.
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(paste0(OUT_STEM, ".svg"), width = w_in, height = h_in)
  draw_flowchart(); invisible(dev.off())
  cat("Saved: ", OUT_STEM, ".svg\n", sep = "")
} else {
  cat("  (SVG skipped: install.packages(\"svglite\") if you need one.)\n")
}

cat("\nCanvas ", CANVAS_W, " x ", CANVAS_H, " mm. The PDF is the one to submit\n",
    "if the journal accepts vector art.\n", sep = "")
