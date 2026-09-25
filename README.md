# Eigenvalue & Eigenvector Calculator

A browser-based **Eigenvalue and Eigenvector Calculator**. Enter any square matrix, and it computes that matrix's eigenvalues, their algebraic and geometric multiplicities, and their corresponding eigenvectors — including complex results and defective (non-diagonalizable) matrices — entirely client-side, with nothing ever sent to a server. When the matrix you enter happens to also be a valid graph Laplacian (weighted or unweighted), the calculator additionally reconstructs and visualizes the graph that matrix represents, as a second, independent feature layered on top of the same input.

**Live Demo:** https://eigenvalueproj.netlify.app/

---

## 1. Project Title

**Eigenvalue & Eigenvector Calculator** (branded in the page header as "CalcuSlay"). It is a static, browser-based web application — there is no installation, account, or backend involved. A user opens the page, types in a matrix, and gets back its full eigen-decomposition.

---

## 2. System Overview

The system solves one core mathematical problem: given a real square matrix `A`, find every eigenvalue `λ` and its associated eigenvector(s) `v` satisfying `A v = λ v`, along with each eigenvalue's algebraic multiplicity (how many times it repeats as a root) and geometric multiplicity (how many independent eigenvectors it actually has).

This is intended for:

- Students and instructors working through linear algebra by hand, who want to check their own eigenvalue/eigenvector calculations.
- Anyone who needs a quick, no-install eigen-decomposition of a matrix — including matrices too large to comfortably solve by hand (up to 100×100).
- Anyone exploring graph theory who wants to see the graph a Laplacian matrix represents, without drawing it by hand.

At a high level, the calculator works in two independent stages that both start from the same matrix the user typed in:

1. **Eigen-decomposition** — always runs, for any square matrix. Produces eigenvalues, multiplicities, and eigenvectors.
2. **Graph analysis** — runs in parallel, independently. If the matrix also happens to satisfy the mathematical definition of a graph Laplacian (see [Section 9](#9-laplacian-detection)), the calculator additionally reconstructs the graph's vertices and edges (with weights, if the Laplacian is weighted) and draws it. If the matrix is not a Laplacian, this second stage simply does not activate — the eigenvalue results are entirely unaffected either way.

The full mathematical implementation behind stage 1 (the QR eigenvalue algorithm, null-space eigenvector computation) and the mathematical portion of stage 2 (Laplacian validation, edge/weight extraction) is documented in detail in **`MATHEMATICAL_CORE.md`**. This document explains the system as a whole — what it does, how the pieces connect, and how to use it — and points to `MATHEMATICAL_CORE.md` wherever the mathematics itself needs full explanation.

---

## 3. Major Features

All of the following are implemented in the current version of the project:

- **Matrix size input** — choose any whole-number size from 1×1 to 100×100.
- **Matrix generation** — clicking Create Matrix (or pressing Enter in the size field) builds an empty input grid of that size.
- **Matrix value input** — type any real number (integers, decimals, negatives) into each cell.
- **Eigenvalue calculation** — for any square matrix, real or complex.
- **Eigenvector calculation** — a basis of eigenvector(s) for every distinct eigenvalue.
- **Algebraic multiplicity** — how many times each eigenvalue repeats as a root.
- **Geometric multiplicity** — the dimension of each eigenvalue's eigenspace.
- **Defective matrix detection** — an eigenvalue is flagged "Defective" when its geometric multiplicity is less than its algebraic multiplicity.
- **Complex eigenvalues** and **complex eigenvectors** — displayed in `a + bi` form.
- **Laplacian detection** — automatic recognition of whether the entered matrix is a valid graph Laplacian.
- **Graph generation** — reconstructing vertices and edges directly from a valid Laplacian.
- **Graph classification** — automatic labeling of recognizable shapes (complete, path, cycle, star, disconnected, empty, single vertex).
- **Unweighted graphs** — the original, still-fully-supported case where every edge has weight 1.
- **Weighted graphs** — automatic detection of non-uniform edge weights, with each weight extracted directly from the Laplacian.
- **Weighted edge visualization** — edge weight is shown both as an on-graph numeric label and, for recognized/small graphs, as the edge's visual length (see [Section 12](#12-weighted-graphs)).
- **Graph interaction** — hover/keyboard-focus highlighting of vertices and edges, click-to-reveal degree (and weighted degree).
- **Graph information** — a live vertex/edge count and graph-type label under the drawing.
- **Copy Results** — copies the input matrix, eigenvalues, multiplicities, and eigenvectors as plain text.
- **Copy Graph Info** — copies vertices, edges, graph type, and (for weighted graphs) each edge's weight as plain text.
- **SVG graph download** — saves the current graph drawing as a standalone `.svg` file.
- **Load Example** — one click loads a ready-made 2×2 matrix.
- **Clear/reset** — empties the current matrix's values without changing its size.
- **Keyboard interaction** — Tab/Enter/Arrow-key navigation between matrix cells, Enter-to-create, Enter-to-calculate.
- **Responsive interface** — usable from a ~375px-wide phone screen up through desktop.
- **Accessibility features** — skip link, ARIA labels/roles, live regions, visible focus states, keyboard-operable graph.
- **Large matrix handling** — an informational notice at 50×50 and above; full support up to 100×100.

---

## 4. Technology Stack

- **HTML5** — page structure (`index.html`).
- **CSS3** — all styling, responsive layout, and light/dark theming (`style.css`), including `prefers-color-scheme` and `prefers-reduced-motion` media queries.
- **JavaScript (vanilla, ES5-style, no frameworks or build step)** — all interactivity, the entire eigenvalue engine, and the entire graph engine (`script.js`).
- **Inline SVG** — the graph visualization is built as SVG DOM elements directly by JavaScript (no charting or graphics library).

There is **no backend, server, or database** — the whole application is static files served to a browser, and every calculation happens client-side. No external mathematics, linear-algebra, or graphing library is used anywhere; the QR eigenvalue algorithm, the null-space solver, and the graph logic are all original code in `script.js` (see `MATHEMATICAL_CORE.md`). The project does **not** run Python, or any other server-side language, in the browser — it is JavaScript from top to bottom.

---

## 5. System Architecture

**Eigenvalue / eigenvector path:**

```
User
  ↓
HTML Interface (index.html: matrix-size field, matrix grid, Calculate button)
  ↓
JavaScript Event Handling (script.js: click/keydown listeners)
  ↓
Matrix Input (readMatrix(): validated n×n array of real numbers)
  ↓
Mathematical Engine (computeEigen → eigenvaluesQR → groupEigenvalues → nullspace)
  ↓
Eigenvalues / Eigenvectors / Multiplicities
  ↓
Results (rendered into the page)
```

**Graph path (runs independently, from the same matrix):**

```
Matrix
  ↓
Laplacian Analysis (laplacianAnalysis: validity + edge/weight extraction)
  ↓
Graph Data (edges, weights, degree, weighted degree, connected components)
  ↓
Graph Classification (classifyGraph: Kn / Pn / Cn / star / disconnected / etc.)
  ↓
Graph Layout (computeLayoutPoints or computeWeightedLayoutPoints: vertex coordinates)
  ↓
SVG Visualization (buildGraphSVG: the actual drawn graph)
  ↓
User Interaction (hover, click, copy, download)
```

**Which parts are which:**

| Layer | Examples | Nature |
|---|---|---|
| **Mathematical** | `eigenvaluesQR`, `nullspace`, `groupEigenvalues`, `laplacianAnalysis` | Pure computation on numbers; no DOM, no pixels. Fully documented in `MATHEMATICAL_CORE.md`. |
| **Data processing** | `computeDegrees`, `computeWeightedDegrees`, `countComponents`, `classifyGraph` | Derives structured graph data (degree, components, shape label) from the mathematical output above; still no pixels. |
| **Presentation / visualization** | `computeLayoutPoints`, `computeWeightedLayoutPoints`, `desiredEdgeLength`, `buildGraphSVG`, all of `style.css`, `index.html` | Decides *where a vertex is drawn* and *how it looks*; consumes the data above but computes nothing new mathematically. |

---

## 6. Complete User Workflow

What actually happens, in order, when someone uses the calculator:

1. The user types a whole number (1–100) into the **Matrix size (n)** field.
2. On **Create Matrix** click (or **Enter** in that field), the size is validated (`createMatrix()`); an invalid size shows an error and creates nothing.
3. If valid, an empty `n × n` grid of input cells is generated (`buildMatrixGrid`).
4. The user types a real number into every cell (Tab / Enter / Arrow keys move between cells).
5. The user clicks **Calculate Eigenvalues & Eigenvectors** (or presses Enter in the last cell, which triggers the same action).
6. `runCalculation()` reads and validates every cell (`readMatrix()`) — any empty or non-numeric cell blocks calculation with an error message.
7. The validated matrix is passed to `computeEigen(matrix, n)`.
8. Internally, `eigenvaluesQR` reduces the matrix to Hessenberg form and runs the shifted-QR algorithm to find every eigenvalue (real or complex).
9. `groupEigenvalues` clusters numerically-equal roots — cluster size = **algebraic multiplicity**.
10. (Algebraic multiplicities are now known, from step 9.)
11. For each distinct eigenvalue, `nullspace(matrix, n, eigenvalue)` computes a basis for `(A − λI)`'s null space — that basis **is** the eigenvector set for that eigenvalue.
12. (Eigenvectors are now known, from step 11.)
13. **Geometric multiplicity** = the number of eigenvectors returned in step 11 for that eigenvalue.
14. Wherever geometric multiplicity is less than algebraic multiplicity, that eigenvalue is flagged **Defective**.
15. `renderResults()` formats and displays the input matrix, every eigenvalue, its multiplicities, and its eigenvector(s).
16. Independently, `renderGraphSection()` passes the *same* matrix to `laplacianAnalysis()`, checking the Laplacian rules (see [Section 9](#9-laplacian-detection)).
17. If valid, `laplacianAnalysis` returns the extracted edges (with weights) and whether the graph is weighted.
18. `classifyGraph()` assigns a structural type label (Kₙ, Pₙ, Cₙ, star, disconnected, etc.), if one applies unambiguously.
19. `computeLayoutPoints` (unweighted) or `computeWeightedLayoutPoints` (weighted) computes where each vertex should be drawn.
20. `buildGraphSVG()` renders the actual SVG: vertices, edges, and (for a weighted graph) weight labels.
21. The user can now hover/click the graph, click **Copy Results** / **Copy Graph Info**, or click **Download Graph** — all working from the same computed results, with no recalculation needed.

If step 16's validation fails at any point, steps 17–20 simply do not run: the Graph Visualization card shows an explanatory message instead, and the eigenvalue results from steps 6–15 are displayed exactly as normal, completely unaffected.

---

## 7. Matrix Input System

- **Supported size:** any whole number from **1×1 to 100×100** (`MIN_N = 1`, `MAX_N = 100` in `script.js`).
- **Create Matrix workflow:** the size field does **not** regenerate the matrix as you type — you must click **Create Matrix**, or press **Enter** while the field is focused. This is intentional: it stops the grid from being rebuilt on every keystroke while typing a larger number like `100`.
- **Numeric input:** integers, decimals (`0.75`), and negative values (`-2.5`) are all accepted in every cell. Validated against the pattern `^-?(\d+\.?\d*|\.\d+)$`.
- **Empty / invalid values:** if any cell is left empty, or contains anything that doesn't parse as a valid real number, clicking Calculate shows an error and focuses the offending cell — calculation does not proceed.
- **Large matrix warning:** creating a matrix of **50×50 or larger** shows the notice *"Large matrix detected. Calculations may take longer depending on your device."* (`LARGE_MATRIX_WARN_THRESHOLD = 50`) — verified: the notice is present at exactly 50×50 and absent at 49×49.
- **Matrix generation behavior:** the grid is built fresh each time `Create Matrix` runs (values are cleared unless explicitly repopulated, as `Load Example` does). Navigating between cells supports Tab, Enter (moves right, then down, then triggers Calculate from the last cell), and the arrow keys.

---

## 8. Eigenvalue and Eigenvector System

This section explains the workflow at a system level; the full algorithmic implementation is in `MATHEMATICAL_CORE.md`.

- **Eigenvalue calculation** uses a shifted-QR algorithm (Hessenberg reduction + implicit-shift QR steps via complex Givens rotations, with deflation) — the same family of method production numerical libraries use, chosen because it handles repeated and clustered eigenvalues far more robustly than expanding a characteristic polynomial and finding its roots.
- **Eigenvalue grouping** clusters numerically-close roots (a repeated eigenvalue emerges from the solver as several close-but-not-identical floating-point values) into one reported value per distinct eigenvalue.
- **Algebraic multiplicity** is the size of that cluster — how many times the eigenvalue repeats as a root.
- **Eigenvector calculation** solves `(A − λI)v = 0` via complex Gaussian elimination, for each distinct eigenvalue, returning a basis for its eigenspace.
- **Geometric multiplicity** is the number of vectors in that basis — the eigenspace's dimension.
- **Repeated eigenvalues:** algebraic multiplicity > 1. These may or may not be defective (see below).
- **Defective matrices:** whenever an eigenvalue's geometric multiplicity is *less than* its algebraic multiplicity, the matrix does not have a full set of independent eigenvectors for that eigenvalue (it is not diagonalizable with respect to that eigenvalue), and the UI marks it **Defective**. Example: `[[2,1],[0,2]]` has eigenvalue `2` with algebraic multiplicity 2 but only one independent eigenvector (geometric multiplicity 1) — a classic Jordan-block example.
- **Complex eigenvalues / eigenvectors:** the solver runs in complex arithmetic throughout, so a matrix like `[[0,-1],[1,0]]` (a 90° rotation) correctly returns the complex conjugate pair `λ = i` and `λ = −i`, each with its own complex eigenvector — no special-casing is needed.

See **`MATHEMATICAL_CORE.md`** (Sections 3–9) for the exact algorithm, tolerances, shift-selection strategy, and the full source code.

---

## 9. Laplacian Detection

A matrix the user enters is treated as a valid graph Laplacian only when **every one** of these checks passes, in this order (`laplacianAnalysis()`):

1. **Every entry is a finite real number.**
2. **The matrix is symmetric** (`L[i][j] = L[j][i]`, within a small tolerance).
3. **Every off-diagonal entry is zero or negative** (a positive off-diagonal entry is never valid).
4. **Every diagonal entry is non-negative.**
5. **Every row sums to zero**, within a numerical tolerance that scales with the matrix size.

The check for "square matrix" is implicit — the UI only ever builds `n × n` grids, so this is always satisfied by construction. The numerical tolerance for every comparison is `1e-6` (`LAPLACIAN_TOL`).

**A. If the matrix IS a valid Laplacian:** the graph is reconstructed automatically and the Graph Visualization card shows the drawing, vertex/edge counts, and (where recognizable) a shape label — see [Section 10](#10-graph-system).

**B. If the matrix is NOT a valid Laplacian:** the Graph Visualization card shows a specific, human-readable explanation of *which* rule failed — for example:

> *"Graph visualization is unavailable: row 2 sums to 3, not 0; every row of a Laplacian matrix must sum to zero."*

> *"Graph visualization is unavailable: the matrix is not symmetric (L[1][2] = -2 but L[2][1] = -1), so it cannot represent an undirected graph."*

**In both cases A and B, the eigenvalue/eigenvector calculation is completely unaffected** — it always runs, on every valid square matrix, whether or not that matrix is a Laplacian. Not every matrix produces a graph; the "Load Example" matrix `[[2,1],[0,3]]`, for instance, is not symmetric and therefore never produces one, while its eigenvalues (`2` and `3`) are still calculated and shown normally.

---

## 10. Graph System

Once `laplacianAnalysis` confirms a matrix is a valid Laplacian, it returns the graph as structured data — never as pixels or coordinates at this stage:

- **Vertices:** row/column index `i` (0-indexed in the matrix) becomes vertex `i + 1`.
- **Edges:** for every off-diagonal entry `L[i][j] < 0` (checked only in the upper triangle, since the matrix is already known symmetric), an edge exists between vertex `i+1` and vertex `j+1`.
- **Weights:** each such edge's weight is `|L[i][j]|`.
- **Degree:** the number of edges touching a vertex (`computeDegrees`).
- **Weighted degree:** the sum of the weights of the edges touching a vertex (`computeWeightedDegrees`) — only computed/shown when the graph is weighted; for an unweighted graph it would equal plain degree anyway.
- **Connected components:** a union-find pass (`countComponents`) over the vertex set, used to detect a disconnected graph.
- **Graph classification:** `classifyGraph` labels the graph Kₙ / Pₙ / Cₙ / star / disconnected / empty / single-vertex when the degree sequence and edge count unambiguously match one of those shapes; otherwise no shape label is invented.

**MATHEMATICAL GRAPH DATA vs. GRAPH VISUALIZATION — the important distinction:**

> Everything above (vertices, edges, weight = `|L[i][j]|`, degree, weighted degree, component count, shape label) is **mathematical graph data** — plain numbers and strings, computed once from the validated matrix.
>
> **None of it is where the graph is drawn on screen.** Turning that data into vertex coordinates, edge lengths, and an actual SVG picture is a separate, later step ([Section 13](#13-graph-layout-system)) that *reads* this data but never changes it. A weighted graph's edge weight is always exactly `|L[i][j]|`, regardless of how long that edge happens to be drawn.

---

## 11. Unweighted Graphs

An unweighted graph is simply the special case where every extracted edge weight equals exactly `1` (`weighted === false`). The calculator recognizes and draws these shapes with a dedicated layout for each, rather than putting every graph on one generic circle:

| Shape | How it's drawn |
|---|---|
| **Path graph** (Pₙ) — e.g. **P₄** | A straight horizontal line, vertices ordered left-to-right by their actual connection order along the path (not necessarily numeric order). |
| **Cycle graph** (Cₙ) — e.g. **C₄** | 4-vertex cycles use an axis-aligned square (one vertex per corner); other sizes use a regular polygon. |
| **Star graph** (K₁,ₙ₋₁) — e.g. **K₁,₃** | The hub (degree `n−1`) sits dead center; every leaf is placed evenly around it. |
| **Complete graph** (Kₙ) — e.g. **K₃**, **K₄** | K₄ uses the classic "triangle with a center vertex" arrangement; other complete graphs use a symmetric regular polygon. |
| **Disconnected graph** | Each connected component gets its own grid cell in the drawing, laid out independently with the same shape rules above — components are never mixed onto one shared shape. |
| **Generic / unrecognized** | A symmetric circular arrangement in vertex-number order, as a reliable fallback for any shape that doesn't match one of the above. |

Path and cycle vertex ordering is derived from the graph's actual connections (a chain-tracing walk), not from raw vertex numbering, so the drawing is always a clean, uncrossed line or ring.

Verified: `K₃` (`[[2,-1,-1],[-1,2,-1],[-1,-1,2]]`) → 3 vertices, 3 edges, **Complete Graph K₃**. `P₄` (`[[1,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,1]]`) → 4 vertices, 3 edges, **Path Graph P₄**. `K₁,₃` → 4 vertices, 3 edges, **Star Graph K₁,₃**. `C₄` → 4 vertices, 4 edges, **Cycle Graph C₄**. `K₄` (`[[3,-1,-1,-1],[-1,3,-1,-1],[-1,-1,3,-1],[-1,-1,-1,3]]`) → 4 vertices, 6 edges, **Complete Graph K₄**, eigenvalues `0, 4, 4, 4`.

---

## 12. Weighted Graphs

A **weighted** graph is detected the moment any single extracted edge weight differs from `1` — there is no separate toggle; the matrix itself determines it (`laplacianAnalysis`'s `weighted` flag).

- **Detection:** while extracting edges, if any weight `w = |L[i][j]|` differs from `1` by more than the numerical tolerance, the whole graph is flagged weighted.
- **Weight extraction:** identical mechanism as the unweighted case — `w_ij = -L_ij`, for every negative off-diagonal entry. An unweighted graph is not a different code path; it is simply the case where every extracted weight is `1`.
- **Mathematical meaning:** the weight is graph data (a distance, cost, or capacity the matrix encodes) — it has nothing to do with eigenvalues. **The calculator never uses eigenvalues to determine anything about the graph or its visualization** — the two pipelines are entirely independent (see [Section 5](#5-system-architecture)).
- **Weighted degree:** the sum of a vertex's incident edge weights, shown alongside plain degree only for weighted graphs (e.g. "Vertex 1 — Degree: 3, Weighted Degree: 4.").
- **How weights are displayed:** every weighted edge gets a numeric label drawn beside it (offset from the line so it doesn't sit on top of it), and hovering an edge shows a tooltip with the exact weight.

**How weight affects the visualization — edge length:**

> **Smaller edge weight → shorter visual edge. Larger edge weight → longer visual edge.**

This mapping (`desiredEdgeLength`) is a bounded min-max normalization against the graph's own weight range — a handful of very large or very small weights cannot blow up or collapse the drawing, and the exact pixel length is not meant to equal the weight number; only its length *relative to the other edges in the same graph* is meaningful. If every edge in a weighted graph happens to share the same weight, there is no relative difference to show, so it falls back to the plain equal-length layout automatically.

Behavior differs by recognized shape (see [Section 13](#13-graph-layout-system) for the exact functions):

- **Weighted paths** — still one straight line; each segment's length along that line scales with its weight.
- **Weighted cycles** — still a ring of the same radius; vertices are spaced *unevenly around the ring* (angular position scales with weight) so lower-weight edges pull their endpoints closer together along the ring and higher-weight edges push them apart — this keeps the ring clean and uncrossed at any size.
- **Weighted stars** — the hub stays centered; each leaf's *distance from the hub* is exactly its edge's weight (normalized), rather than every leaf sitting on one fixed circle.
- **Small weighted complete/generic graphs** (≤ 14 vertices in one connected component) — positioned with a lightweight, bounded, damped force simulation: edges act as springs pulled toward a weight-proportional rest length, vertices repel each other so nothing overlaps, seeded from the same shape-aware starting layout an unweighted graph would use.
- **Large weighted graphs** (> 14 vertices, complete or unrecognized shape) — kept at the plain equal-length layout; once a graph is already dense enough to be visually hard to trace, stretching edges by weight would add noise, not clarity. Weight labels and hover tooltips remain the reliable way to read exact weights at that size.

Verified: the weighted-K₄ matrix below produces edges 1–2=2, 1–3=1, 1–4=4, 2–3=3, 2–4=3, 3–4=1, and the rendered/downloaded SVG's measured pixel lengths were strictly ordered by weight (shortest for the two weight-1 edges, longest for the weight-4 edge) — see [Section 21](#21-important-verified-mathematical-examples).

---

## 13. Graph Layout System

Layout is entirely separate from graph *data* (Section 10) — it only decides pixel coordinates, never which edges exist.

**Unweighted path:** `computeLayoutPoints(n, edges, size, nodeRadius)` — dispatches each connected component (`componentsList`) to `detectComponentShape` (path / cycle / star / complete / generic) and then `layoutComponent`, which places vertices per the table in [Section 11](#11-unweighted-graphs). Multiple components are each given their own cell in a grid.

**Weighted path:** `computeWeightedLayoutPoints(n, edges, size, nodeRadius)` — the weighted-graph counterpart, following the same component/grid structure, but with per-shape logic:

| Component shape | Weighted layout function |
|---|---|
| Path | `layoutWeightedPath` — analytic placement along a line; segment length ∝ weight. |
| Cycle | `layoutWeightedCycle` — analytic placement around a fixed-radius ring; angular spacing ∝ weight. |
| Star | `layoutWeightedStar` — hub centered; leaf radius from hub ∝ weight. |
| Complete / generic, ≤ 14 vertices | `layoutComponent` (shape-aware seed) followed by `forceDirectedLayout` (bounded, damped spring simulation; rest length per edge from `desiredEdgeLength`). |
| Complete / generic, > 14 vertices | `layoutComponent` only (plain equal-length layout; no weight-based distortion). |
| Any component, uniform weight | The whole graph defers entirely to `computeLayoutPoints` (the unweighted layout) — there is no relative weight difference to visualize. |

`desiredEdgeLength(weight, minWeight, maxWeight, minLen, maxLen)` performs the actual weight → pixel-length mapping: min-max normalize the weight into `[0, 1]`, then map that into a `[minLen, maxLen]` pixel range — bounded regardless of how extreme the input weights are.

**Component handling** (both paths): a graph is split into connected components (`componentsList`/`countComponents`); a single component fills the whole drawing area, while multiple components are each assigned their own cell in an automatically-sized grid, so disconnected graphs never draw one component's edges crossing into another's space.

**Unweighted layouts are entirely unaffected by any of this** — `buildGraphSVG` only ever calls `computeWeightedLayoutPoints` when `weighted === true`; an unweighted graph's layout call, and therefore its drawing, is identical to before weighted-graph support existed.

---

## 14. Graph Visualization and Interaction

- **SVG generation** (`buildGraphSVG`) builds the drawing as real SVG DOM elements (`<svg>`, `<line>`, `<circle>`, `<text>`, `<g>`) — not a canvas bitmap or an external image.
- **Vertices** are circles, numbered (for small/medium graphs — see below), colored from the page's own CSS custom properties so the graph matches the current light/dark theme.
- **Edges** are straight lines between the vertex coordinates from Section 13.
- **Labels:** vertex numbers are shown permanently for graphs up to ~45 vertices; above that, only on hover/focus, to avoid clutter.
- **Weight labels:** for a weighted graph, each edge's weight is drawn as a small boxed number offset to one side of the edge — under the same ~45-vertex threshold as vertex labels; above it, the weight is available via hover only (see [Section 20](#20-testing-and-verified-behavior) for the verified density note text).
- **Tooltips:** every edge and vertex has a native `<title>` tooltip (e.g. *"Edge 1-2 / Weight: 2"*, *"Vertex 1 (degree 3)"*), shown on hover by the browser.
- **Hover behavior:** hovering a vertex highlights it and every edge touching it; hovering an edge (or its weight label) highlights that edge and its label together.
- **Focus behavior:** Tab-focusing a vertex produces the same highlight as hovering it, for keyboard users.
- **Click behavior:** clicking a vertex (or pressing Enter/Space while it's focused) reveals its degree in a line below the graph.
- **Degree display:** for an unweighted graph, "Vertex 1 — degree 3." For a **weighted** graph: "Vertex 1 — Degree: 3, Weighted Degree: 4." — both numbers shown together, since they answer different questions.
- **Copy Graph Info:** copies vertex count, edge count, graph type, and the full edge list to the clipboard — including each edge's weight, for a weighted graph (e.g. `1-2: weight 2`).
- **Download Graph SVG:** serializes the currently-displayed SVG (exactly as rendered, including weighted layout and weight labels) into a standalone `laplacian-graph.svg` file.

None of these interactions are required to read the graph — vertex/edge counts, the graph-type label, and the drawing itself (including weight labels, where shown) are always visible without hovering or clicking anything.

---

## 15. User Interface

The page (`index.html`) is a single static document with these sections, top to bottom:

- **Header** — hero card with the "CalcuSlay" brand mark and the page title/subtitle.
- **"How it works"** — a 4-step at-a-glance summary card.
- **Matrix Size section** — the size field and **Create Matrix** button.
- **Matrix Entries section** — the generated input grid, plus **Calculate Eigenvalues & Eigenvectors**, **Load Example**, and **Clear Matrix** buttons, and a "Calculating…" indicator.
- **Results section** (hidden until a calculation completes):
  - **Input Matrix** card — a static, formatted redisplay of what was entered.
  - **Eigenvalues & Eigenvectors** card — one card per distinct eigenvalue, with a **Copy Results** button.
  - **Graph Visualization** card — the SVG drawing (or an explanatory unavailable-message), a "Weighted Graph" badge when applicable, **Copy Graph Info** and **Download Graph** buttons, and a density note for large/dense graphs.
- **Feedback/messages** — a single error box (`role="alert"`) for input problems, and a large-matrix notice for 50×50+.
- **Footer** — a one-line reminder that all calculations run locally.

This document does not reproduce `style.css` rule-by-rule — see that file directly for exact colors, spacing, and typography. What matters structurally is that every interactive element above has a stable `id` that `script.js` binds to directly (no framework, no data-binding layer).

---

## 16. Responsive and Accessibility Behavior

- **Desktop layout:** full-width cards in a centered container.
- **Mobile layout:** CSS breakpoints at `700px` and `560px` reflow the matrix grid, button rows, and graph card for narrow screens; verified down to ~375px wide with no horizontal overflow and no label clipping.
- **Keyboard navigation:** Tab/Shift+Tab moves through controls; inside the matrix grid, Tab, Enter, and the Arrow keys all move between cells (`onCellKeydown`).
- **Enter key:** pressing Enter in the size field triggers Create Matrix; pressing Enter in the last matrix cell triggers Calculate — both verified.
- **Focus behavior:** all interactive elements have visible `:focus-visible` outlines; the graph's vertices are individually focusable (`tabindex="0"`) and keyboard-activatable (Enter/Space).
- **ARIA / accessible labels:** a skip-link to `#main-content`; `aria-labelledby` on major sections; `role="alert"` on the error box; `role="status"` + `aria-live="polite"` on the calculating indicator and large-matrix notice; `role="img"` + a descriptive `aria-label` on the graph container (which updates to mention weights when the graph is weighted); `aria-live="polite"` on the vertex-degree readout; each matrix cell has its own `aria-label` ("Row 2, column 3"); each SVG vertex/edge has a native `<title>` tooltip as well as `role="button"` + `aria-label` on vertices.
- **Reduced motion:** a `prefers-reduced-motion: reduce` media query is present in `style.css` for users who request less animation.
- **Theme:** light/dark theme follows the OS-level `prefers-color-scheme` automatically — including the graph's own colors, which are read from CSS custom properties at render time so a graph rendered in dark mode is dark, not just the page around it.

Only behavior actually present in the current `index.html`/`style.css`/`script.js` is listed above.

---

## 17. Project Structure

```text
Eigenvalue & Eigenvector Calculator/
├── .gitattributes
├── index.html
├── script.js
├── style.css
├── README.md
└── MATHEMATICAL_CORE.md
```

| File | Purpose / responsibility | Does NOT control |
|---|---|---|
| **`index.html`** | Page structure: every section, button, input, and container the calculator uses; the fixed pieces of markup `script.js` fills in and reads from. | Styling, colors, layout math, or any calculation. |
| **`style.css`** | All visual styling: colors, spacing, typography, responsive breakpoints, light/dark theming, focus states, reduced-motion behavior. | Any calculation, any DOM structure, any graph data. |
| **`script.js`** | Everything interactive: matrix grid generation, input validation, the entire eigenvalue engine, the entire graph engine (validation, data, layout, SVG rendering), and all button/keyboard behavior. | Visual appearance (delegates to CSS classes it applies, but doesn't define how those classes look). |
| **`README.md`** *(this file)* | Complete system-level documentation: what the app is, how it works end to end, how to use and deploy it. | Line-by-line mathematical derivations — see `MATHEMATICAL_CORE.md`. |
| **`MATHEMATICAL_CORE.md`** | The mathematical/algorithmic reference: full source of the eigenvalue engine and the graph-data mathematics, with verified worked examples. | UI structure, styling, deployment, or anything not directly mathematical. |
| **`.gitattributes`** | Line-ending normalization for the repository. | Anything runtime-related. |

---

## 18. Data Flow

**Mathematical calculation:**

```
Matrix input (validated n×n real array, from readMatrix())
  → computeEigen(matrix, n)
      → eigenvaluesQR(matrix, n)        [Hessenberg reduction + shifted QR]
      → groupEigenvalues(roots)         [clusters -> algebraic multiplicity]
      → nullspace(matrix, n, λ)         [per distinct eigenvalue -> eigenvectors]
  → { value, multiplicity, vectors }  per eigenvalue
  → formatted results rendered on the page
```

**Graph generation:**

```
Matrix (the same array)
  → laplacianAnalysis(matrix, n)
      → { ok:false, reason }                         (not a Laplacian -- no graph)
      → { ok:true, edges:[[i,j,weight],...], weighted }
  → computeDegrees(n, edges)  /  computeWeightedDegrees(n, edges)
  → countComponents(n, edges)
  → classifyGraph(n, edges, degrees)
  → computeLayoutPoints(...)  or  computeWeightedLayoutPoints(...)
  → buildGraphSVG(...)
  → user interaction (hover / click / copy / download)
```

---

## 19. Important Functions

A high-level reference; full implementations (for the mathematical functions) are in `MATHEMATICAL_CORE.md`.

| Function | File | Purpose | Category |
|---|---|---|---|
| `eigenvaluesQR` | `script.js` | Shifted-QR eigenvalue solver (Hessenberg + Givens rotations, deflation). | Mathematical |
| `groupEigenvalues` | `script.js` | Clusters numerically-equal roots; determines algebraic multiplicity. | Mathematical |
| `nullspace` | `script.js` | Solves `(A − λI)v = 0`; returns eigenvectors, whose count is geometric multiplicity. | Mathematical |
| `computeEigen` | `script.js` | Public entry point tying the three functions above together. | Mathematical |
| `laplacianAnalysis` | `script.js` | Validates a matrix as a (weighted or unweighted) Laplacian; extracts edges and weights. | Mathematical / Data |
| `computeDegrees` | `script.js` | Plain vertex degree from the edge list. | Data processing |
| `computeWeightedDegrees` | `script.js` | Weighted vertex degree from the edge list. | Data processing |
| `countComponents` | `script.js` | Connected-component count (union-find). | Data processing |
| `classifyGraph` | `script.js` | Structural shape classification (Kₙ, Pₙ, Cₙ, star, disconnected, etc.). | Data processing |
| `computeLayoutPoints` | `script.js` | Unweighted, shape-aware vertex coordinate layout. | Presentation |
| `computeWeightedLayoutPoints` | `script.js` | Weighted vertex coordinate layout (analytic + bounded force simulation). | Presentation |
| `desiredEdgeLength` | `script.js` | Normalizes a weight into a bounded pixel-length range. | Presentation |
| `buildGraphSVG` | `script.js` | Builds the actual SVG graph drawing and its interaction handlers. | Presentation |
| `runCalculation` | `script.js` | UI workflow: reads the matrix, calls `computeEigen` and the graph pipeline, renders results. | UI workflow |
| `readMatrix` / `buildMatrixGrid` | `script.js` | Reads/validates cell input; generates the input grid. | UI workflow |
| `copyResults` / `copyGraphInfo` / `downloadGraph` | `script.js` | Clipboard and file-export actions. | UI workflow |

---

## 20. Testing and Verified Behavior

Every result below was obtained by actually running the current code (either the exact extracted mathematical functions in Node.js, or the live page in a headless browser) — none of it is invented or assumed.

| Test | Verified result |
|---|---|
| 1×1 (`[7]`) | `λ = 7`, algebraic/geometric multiplicity 1. |
| Distinct eigenvalues (`[[2,1],[0,3]]`, the Load Example matrix) | `λ = 2`, `λ = 3`, each multiplicity 1. |
| Repeated, non-defective (`[[2,0],[0,2]]`) | `λ = 2`, algebraic 2, geometric 2 — not defective. |
| Defective (`[[2,1],[0,2]]`) | `λ = 2`, algebraic 2, geometric **1** — flagged Defective. |
| Complex eigenvalues (`[[0,-1],[1,0]]`) | `λ = −i` and `λ = i`, each multiplicity 1. |
| Identity matrix (`3×3`) | `λ = 1`, algebraic/geometric multiplicity 3. |
| Zero matrix (`3×3`) | `λ = 0`, algebraic/geometric multiplicity 3. |
| Diagonal matrix (`[[5,0,0],[0,-2,0],[0,0,7]]`) | `λ = −2, 5, 7`, each multiplicity 1. |
| K₃ Laplacian | `λ = 0, 3, 3`; vertices=3, edges=3, **Complete Graph K₃**. |
| P₄ Laplacian | `λ = 0, 0.585786, 2, 3.414214`; vertices=4, edges=3, **Path Graph P₄**. |
| K₁,₃ Laplacian | `λ = 0, 1, 1, 4`; vertices=4, edges=3, **Star Graph K₁,₃**. |
| C₄ Laplacian | `λ = 0, 2, 2, 4`; vertices=4, edges=4, **Cycle Graph C₄**. |
| K₄ Laplacian | `λ = 0, 4, 4, 4`; vertices=4, edges=6, **Complete Graph K₄**. |
| **Weighted K₄** | `λ = 0, 5.763932, 10.236068, 12`; vertices=4, edges=6, **Weighted Complete Graph K₄**; edge lengths in the rendered/downloaded SVG strictly ordered by weight. |
| **Weighted P₄** | `λ = 0, 1, 3, 8`; vertices=4, edges=3, **Weighted Path Graph P₄**; still a straight line, segment lengths ordered by weight. |
| **Weighted C₄** | `λ = 0, 2.855148, 6.718356, 10.426496`; vertices=4, edges=4, **Weighted Cycle Graph C₄**. |
| **Weighted star** | `λ = 0, 1.202521, 2.612831, 10.184648`; vertices=4, edges=3, **Weighted Star Graph K₁,₃**; leaf distance from hub ordered by weight. |
| **Large weighted cycle** (50 vertices, weights cycling 1–6) | Renders as a clean, uncrossed ring with uneven vertex spacing reflecting weight; correctly classified **Weighted Cycle Graph C₅₀**; permanent weight labels correctly suppressed above ~45 vertices, with a density note and working hover tooltip. |
| Disconnected graph (two separate edges) | `λ = 0, 0, 2, 2`; vertices=4, edges=2, **Disconnected Graph (2 components)**. |
| Invalid Laplacian (`[[3,-2,-1],[-2,5,0],[-1,0,1]]`) | Eigenvalues still calculated and shown normally; graph rejected with *"row 2 sums to 3, not 0; every row of a Laplacian matrix must sum to zero."* |
| Asymmetric matrix (`[[2,-2,0],[-1,3,-2],[0,-2,2]]`) | Rejected: *"the matrix is not symmetric (L[1][2] = -2 but L[2][1] = -1)..."* |
| Invalid matrix sizes (`0`, `101`) | Both rejected with *"Matrix size must be a whole number between 1 and 100."*; no matrix is created. |
| 100×100 creation | Creates exactly 10,000 input cells; a full calculation completed in a live browser run with zero console errors, producing correctly-formatted results. |
| Enter-to-create | Pressing Enter in the size field builds the matrix grid, same as clicking Create Matrix. |
| Enter-to-calculate | Pressing Enter in the last matrix cell triggers Calculate; verified against a K₃ Laplacian, which produced the correct graph. |
| Graph download | The downloaded `.svg`'s edge coordinates matched the on-screen layout exactly, including weighted edge lengths and weight-label elements. |
| Copy Graph Info | Clipboard content verified to include vertices, edges, graph type, and (for a weighted graph) each edge's weight, in the documented format. |
| Console-error checks | Zero console or page errors across every test above, including the 100×100 run and the 50-vertex weighted graph. |

---

## 21. Important Verified Mathematical Examples

**Weighted K₄** — the canonical example for this project:

```
 7  -2  -1  -4
-2   8  -3  -3
-1  -3   5  -1
-4  -3  -1   8
```

**Verified eigenvalues:**

```
0
5.763932
10.236068
12
```

> `0, 4, 8, 16` is **incorrect** for this matrix — it does not satisfy `det(A − λI) = 0` for `λ = 4`, `8`, or `16`. The set above is the correct, independently-verified result. See `MATHEMATICAL_CORE.md`, Section 14, for the verification method.

Edges: 1–2 = 2, 1–3 = 1, 1–4 = 4, 2–3 = 3, 2–4 = 3, 3–4 = 1 — topology K₄ (6 edges), but drawn with visibly different edge lengths per weight rather than the perfect symmetric triangle-plus-center an unweighted K₄ would use.

**Weighted P₄:**

```
 2  -2   0   0
-2   5  -3   0
 0  -3   4  -1
 0   0  -1   1
```

Edges: 1–2 = 2, 2–3 = 3, 3–4 = 1. Verified eigenvalues: `0, 1, 3, 8`. Drawn as a straight line, 1—2—3—4, with the weight-1 segment shortest and the weight-3 segment longest.

For the full set of additional verified examples (K₃, K₁,₃, C₄, K₄, weighted C₄, weighted star, identity, zero matrix, defective/complex-eigenvalue cases, etc.), see `MATHEMATICAL_CORE.md`, Section 14, and [Section 20](#20-testing-and-verified-behavior) above.

---

## 22. Limitations

- **Numerical, not symbolic:** this calculator uses floating-point numerical methods (shifted QR, Gaussian elimination), not exact symbolic algebra. Results are typically accurate to several decimal places but may show very small rounding artifacts.
- **Clustered eigenvalues:** eigenvalues that are very close together numerically can be harder for any numerical solver — including this one — to separate with full precision.
- **Real-valued matrix input only:** every matrix entry must be a real number; the interface does not accept complex entries directly (only complex *results* are supported).
- **Performance at large sizes:** a 100×100 matrix's calculation time varies with the specific matrix (how quickly its eigenvalues converge under shifted QR) — measured between roughly 2 and 9 seconds across different random 100×100 matrices on the machine used for testing; performance depends on the user's own device and browser.
- **Square matrices only:** the interface only ever builds square grids (1×1 to 100×100); there is no non-square matrix support.
- **Simple graphs only:** graph visualization supports simple, undirected graphs (no self-loops, no multiple edges between the same pair of vertices) — both weighted and unweighted. Directed graphs and multigraphs are not representable as a Laplacian under the rules in Section 9 and will not produce a graph.
- **Edge weights must be positive:** a weight is always `|L[i][j]|` for a negative off-diagonal entry; the calculator never infers a weighted graph from a matrix that fails Laplacian validation.
- **Graph classification is deliberately conservative:** only the specific, unambiguous shapes in Section 11 are labeled. A valid graph that doesn't match one of them is still drawn correctly (with weights, if applicable) but without an invented shape label — an unweighted graph shows only vertex/edge counts in that case, while a weighted graph still shows the generic "Weighted Graph" label, since the weights remain informative on their own.
- **Layout is deterministic but not draggable:** vertices cannot be manually repositioned by the user. Unweighted graphs and the analytic weighted path/cycle/star layouts are fully deterministic (the same matrix always lays out the same way); the bounded force simulation used for small weighted complete/generic graphs is also deterministic (no randomness) but is an approximate physical relaxation, not an exact formula.
- **Large, dense graphs are visually approximate by design:** past a certain size, individual edges become hard to trace by eye even though every one is drawn correctly, permanent labels are intentionally suppressed, and (for weighted graphs) edge length stops reflecting weight in favor of the plain equal-length layout — vertex/edge counts, the graph type, and hover tooltips remain accurate regardless of density.

---

## 23. Deployment

This is a fully static site (HTML, CSS, and JavaScript only) with **no build step, no backend, and no database** — it runs entirely in the visitor's browser. It is currently deployed on **Netlify** at https://eigenvalueproj.netlify.app/.

Required files for deployment: `index.html`, `script.js`, `style.css` (and, if you want the documentation to travel with the deployed copy, `README.md`/`MATHEMATICAL_CORE.md`, though those are not fetched by the page itself). To deploy your own copy on Netlify: push the repository to GitHub, create a new Netlify site connected to it, leave the build command empty, and set the publish directory to the folder containing `index.html` — then deploy. The same approach works on any other static host (GitHub Pages, Vercel, Cloudflare Pages, or simply opening `index.html` directly in a browser / serving it with `python -m http.server`), since there are no server-side requirements anywhere in the project.

---

## 24. Relationship Between Documentation Files

This project intentionally has exactly two Markdown files, each with a distinct, non-overlapping purpose:

- **`README.md`** *(this file)* — **complete documentation of the whole system**: what the application is, every feature, the full user workflow, the UI structure, deployment, and verified test behavior. It explains the mathematics and the graph engine *at a system level* — what each stage does and how it connects to the rest of the app — without reproducing the underlying algorithms line by line.
- **`MATHEMATICAL_CORE.md`** — **mathematical and algorithmic documentation**, purpose-built to also serve as a research/thesis-suitable technical reference. It contains the actual current source code of the eigenvalue engine and the graph-data mathematics, explained in full algorithmic detail, plus a complete, copy-pasteable source-code appendix. It is not a second README and does not cover UI, deployment, or feature lists.
- **The source code itself** (`index.html`, `script.js`, `style.css`) is the ultimate authority. Both documents describe the *current* implementation as of the same revision; if any statement in either document ever conflicts with what the code actually does, the code is correct and the documentation should be treated as needing an update, not the other way around.
