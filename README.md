# Eigenvalue & Eigenvector Calculator

A browser-based calculator that computes the eigenvalues and eigenvectors of any square matrix. Enter a matrix size, fill in the values, and get the eigenvalues, their multiplicities, and their corresponding eigenvectors — all calculated instantly in your browser, with no data ever sent to a server. When the matrix you enter is a graph Laplacian — weighted or unweighted — the calculator also reconstructs and displays the underlying graph automatically, including each edge's weight when the graph is weighted.

**Live Demo:** https://eigenvalueproj.netlify.app/

---

## Overview

This tool lets you:

- Choose the size of a square matrix, from 1×1 up to 100×100
- Enter the matrix values by hand
- Calculate the matrix's eigenvalues
- Calculate the eigenvector(s) associated with each eigenvalue
- See the algebraic and geometric multiplicity of each eigenvalue
- Copy the full set of results as plain text
- Automatically visualize the graph behind a Laplacian matrix, when the matrix you entered is one — including edge weights, when the matrix is a weighted Laplacian

In linear algebra, an eigenvalue and eigenvector pair describes a direction that a matrix only stretches or shrinks, without changing — a property used throughout engineering, physics, computer graphics, and data science (for example, in vibration analysis, stability analysis, and principal component analysis). This calculator handles that computation for you, including matrices whose eigenvalues turn out to be complex numbers.

A graph Laplacian is a special square matrix (`L = D - A`, degree matrix minus adjacency matrix) built from a graph's vertices and edges; its eigenvalues reveal structural properties of the graph (for example, the number of zero eigenvalues equals the number of connected components). This calculator recognizes valid graph Laplacians automatically — both unweighted and weighted — and draws the graph they represent alongside the usual eigenvalue results.

---

## Features

- Square matrix input for any size from 1×1 to 100×100
- Eigenvalue calculation, including complex (non-real) eigenvalues
- Eigenvector calculation for every eigenvalue
- Algebraic multiplicity (how many times an eigenvalue repeats as a root)
- Geometric multiplicity (how many independent eigenvectors an eigenvalue actually has)
- A "Defective" indicator when an eigenvalue's algebraic and geometric multiplicities differ
- Automatic Laplacian detection and graph reconstruction, with a clean, interactive SVG visualization
- Automatic detection of **weighted** graph Laplacians, with each edge's weight displayed directly on the graph
- Automatic graph type labeling for recognizable shapes (complete, path, cycle, star, disconnected, or empty graphs), prefixed with "Weighted" when the graph is weighted
- Hover and keyboard-focus highlighting of graph vertices and their edges, and click-to-reveal vertex degree (and weighted degree, for weighted graphs)
- Downloadable graph image (SVG, including edge weight labels) and one-click copying of the graph's vertex/edge summary (including weights)
- A built-in example matrix you can load with one click
- One-click clearing of the current matrix
- One-click copying of the full eigenvalue/eigenvector results as plain text
- Input validation with clear, specific error messages
- A warning notice for large matrices (50×50 and above)
- Keyboard navigation between matrix cells (Tab, Enter, and arrow keys)
- Responsive layout that works on desktop, tablet, and mobile
- Automatic light/dark theme based on your system preference
- All calculations run locally in your browser — no data is ever uploaded

---

## How to Use

### 1. Choose the Matrix Size

In the **Matrix Size** section, enter a whole number for **n** into the **Matrix size (n)** field. This determines the dimensions of the square matrix you'll fill in.

Examples:

- `2` → a 2×2 matrix
- `3` → a 3×3 matrix
- `4` → a 4×4 matrix

Accepted values are whole numbers from **1 to 100**. Once you've typed a size, click **Create Matrix** (or press **Enter** while the field is focused) to build the matrix grid. The grid is not generated automatically while you type — this is intentional, so that typing a larger number like `100` doesn't rebuild the matrix on every keystroke.

If the value you entered is missing, not a whole number, or outside the 1–100 range, an error message appears and no matrix is created.

### 2. Enter the Matrix Values

Once the matrix grid appears, click into each cell and type its value. Negative numbers and decimals are supported (for example `-2.5` or `0.75`).

For a 2×2 matrix, you might enter:

```
[ 2  1 ]
[ 0  3 ]
```

You can move between cells with **Tab**, **Enter**, or the **arrow keys**. Pressing **Enter** in the last cell of the grid automatically triggers the Calculate button.

### 3. Calculate

Click **Calculate Eigenvalues & Eigenvectors** once every cell has a value. The calculator reads the matrix, runs the calculation, and displays the results below. Every cell must contain a valid number — if any cell is empty or contains an invalid value, an error message appears and the calculator will not run.

### 4. Read the Results

For each distinct eigenvalue found, the results section shows:

- **Eigenvalue (λ)** — the numeric value itself
- **Algebraic Multiplicity** — how many times this eigenvalue is repeated as a root of the matrix
- **Geometric Multiplicity** — how many independent eigenvectors exist for this eigenvalue
- **Eigenvector(s)** — the vector, or basis of vectors, associated with the eigenvalue

If an eigenvalue's geometric multiplicity is lower than its algebraic multiplicity, it is labeled **Defective**, meaning the matrix doesn't have a full set of independent eigenvectors for that eigenvalue.

### 5. Load Example

Click **Load Example** to instantly load a ready-made 2×2 matrix:

```
[ 2  1 ]
[ 0  3 ]
```

This is a quick way to see how the calculator works without typing your own values.

### 6. Clear Matrix

Click **Clear Matrix** to empty every cell in the current matrix grid. The matrix size stays the same — only the values are removed. This also clears any results currently on screen.

### 7. Copy Results

After a calculation completes, click **Copy Results** to copy the input matrix, the eigenvalues (with their multiplicities), and all eigenvectors to your clipboard as plain text, ready to paste elsewhere.

### 8. Graph Visualization (for Laplacian matrices)

If the matrix you calculated is a valid graph Laplacian — unweighted or weighted — a **Graph Visualization** card automatically appears below the eigenvalue results, showing the graph the matrix represents, with edge weights drawn on the graph when it's weighted. If it isn't, the card instead shows a short message explaining specifically why graph visualization isn't available for that matrix — the eigenvalue and eigenvector results above are unaffected either way. See [Laplacian Matrix and Graph Visualization](#laplacian-matrix-and-graph-visualization) below for full details.

---

## Understanding the Results

### Eigenvalues

An eigenvalue λ is a number that satisfies the equation:

```
A v = λ v
```

where **A** is the matrix, **v** is a nonzero vector (the eigenvector), and **λ** is the eigenvalue. In other words, when the matrix is applied to its eigenvector, the result is the same vector scaled by λ — its direction doesn't change.

### Eigenvectors

An eigenvector is the vector **v** in the equation above. It represents a direction that the matrix only stretches, shrinks, or flips — it does not rotate it into a different direction.

### Algebraic Multiplicity

The algebraic multiplicity of an eigenvalue is how many times it appears as a repeated root when solving for the matrix's eigenvalues. For example, an eigenvalue with algebraic multiplicity 2 means the matrix has that same eigenvalue counted twice.

### Geometric Multiplicity

The geometric multiplicity of an eigenvalue is the number of independent eigenvectors associated with it. This can be equal to or less than the algebraic multiplicity, but never greater. When it's less, the calculator labels that eigenvalue as **Defective**.

### Complex Eigenvalues

Not every matrix has purely real eigenvalues — some produce complex results. The calculator fully supports this and displays complex eigenvalues and eigenvectors in the form:

```
a + bi
```

For example, a result might appear as `i` or `-i`, or as a combined value such as `1 + 2i`.

---

## Example

Given the matrix:

```
A =
[ 2  1 ]
[ 0  3 ]
```

The calculator returns two eigenvalues: **2** and **3**, each with algebraic and geometric multiplicity 1, along with their corresponding eigenvectors. This is the same matrix used by the **Load Example** button.

---

## Laplacian Matrix and Graph Visualization

### What Is a Graph Laplacian?

For a simple, undirected, unweighted graph, the **Laplacian matrix** is defined as:

```
L = D - A
```

where **D** is the diagonal degree matrix (each diagonal entry is the number of edges touching that vertex) and **A** is the adjacency matrix (a 1 in position `[i][j]` means an edge connects vertex i and vertex j, 0 otherwise). The result is a symmetric matrix where every row sums to zero. The Laplacian's eigenvalues describe structural properties of the graph — for example, the eigenvalue 0 always appears, and it appears once for every connected component in the graph.

Whenever you calculate a matrix that fits this pattern, the calculator reconstructs the graph directly from the matrix entries and displays it automatically — you don't need to do anything extra. This same idea extends naturally to **weighted** graphs, covered in its own section below.

### How the Calculator Recognizes a Laplacian

A calculated matrix is treated as a valid graph Laplacian only when **all** of the following hold:

- It is square (guaranteed by the matrix size you chose).
- Every entry is a finite real number.
- It is symmetric (`L[i][j] = L[j][i]` for every pair).
- Every off-diagonal entry is zero or negative (a positive off-diagonal entry is never valid).
- Every diagonal entry is non-negative.
- Every row sums to zero (within a small numerical tolerance).

If any of these checks fail, the matrix is not treated as a Laplacian. The eigenvalue and eigenvector results are calculated and shown exactly as normal, and the Graph Visualization card explains specifically why the matrix was rejected, for example:

> "Graph visualization is unavailable: row 2 sums to 3, not 0; every row of a Laplacian matrix must sum to zero."

or:

> "Graph visualization is unavailable: the matrix is not symmetric (L[1][2] = -2 but L[2][1] = -1), so it cannot represent an undirected graph."

These are informational messages, not errors, and the calculator never tries to silently "fix" or reinterpret an invalid matrix — the rest of the calculator continues to work normally either way.

A matrix that passes these checks is automatically classified as either an **unweighted** or a **weighted** Laplacian (see below) — there is no separate toggle to set; the matrix itself determines which one you get.

### How the Graph Is Reconstructed

The graph is built directly from the matrix, never from the eigenvalues:

- Each row/column index becomes a vertex, numbered starting from **1** (row/column 1 in the matrix is vertex 1, and so on).
- For every off-diagonal entry `L[i][j] < 0` (with i < j), one edge is drawn between vertex i and vertex j, and its weight is the entry's absolute value. Because the matrix is symmetric, the mirrored entry `L[j][i]` is not read again, so each edge is created exactly once.
- An off-diagonal entry of exactly `0` means there is no edge between that pair of vertices.

For example, the Laplacian

```
 3 -1 -1 -1
-1  3 -1 -1
-1 -1  3 -1
-1 -1 -1  3
```

reconstructs into 4 vertices, each connected to every other vertex — 6 edges in total (the complete graph K₄), each with weight 1.

### Weighted Graphs

A weighted graph assigns a numerical weight to each edge — for example, a distance, a cost, or a capacity. The calculator recognizes this directly from the matrix, with no separate setting to turn on: whenever the off-diagonal entries of a valid Laplacian aren't all exactly `-1`, the graph is treated as weighted.

For a weighted Laplacian:

- A negative off-diagonal entry represents an edge, and **its absolute value is the edge's weight**. For example, `L[1][2] = -2` means an edge between vertex 1 and vertex 2 with weight `2`.
- An off-diagonal entry of `0` means there is no edge, exactly as in the unweighted case.
- Each diagonal entry equals the sum of the weights of every edge touching that vertex (its **weighted degree**), not just the number of edges.
- The matrix must still be symmetric, and every row must still sum to zero, exactly as with an unweighted Laplacian — these two rules are what make a weighted Laplacian mathematically valid.

For example, the matrix

```
 4 -2 -1 -1
-2  5 -3  0
-1 -3  4  0
-1  0  0  1
```

is a valid weighted Laplacian (every row sums to zero, and it's symmetric) representing 4 vertices and 4 weighted edges: 1–2 with weight 2, 1–3 with weight 1, 1–4 with weight 1, and 2–3 with weight 3. There is no edge between 2–4 or 3–4, because those off-diagonal entries are `0`.

When a weighted Laplacian is detected, the Graph Visualization card:

- Labels the graph as **Weighted Graph**, or, when the underlying shape is confidently recognized, as **Weighted Complete Graph**, **Weighted Path Graph**, **Weighted Cycle Graph**, and so on — the same shape recognition used for unweighted graphs, just with a "Weighted" prefix.
- Draws each edge's weight directly beside it, offset slightly from the line so the number stays readable rather than sitting on top of it.
- **Draws each edge's length to reflect its weight** — a smaller weight draws as a visually shorter edge and a larger weight as a visually longer one, so the relative size of every weight is visible at a glance, not just readable from the labels. This is the one respect in which a weighted graph's layout differs from the unweighted case: an unweighted graph (every edge implicitly weight 1) still uses the fixed, equal-length shape templates described below, while a weighted graph stretches or compresses edges by weight while keeping the same vertices connected. The exact pixel length isn't meant to equal the weight number precisely — only its length *relative to the other edges* is meaningful, normalized so a handful of very large or very small weights can't blow up or collapse the drawing. If every edge happens to share the same weight, there's no relative difference to show, so the graph falls back to the plain equal-length layout, exactly like an unweighted graph.
- Never changes which vertices are connected — the weight only affects how long an edge is drawn, never whether it exists.

#### Degree vs. Weighted Degree

These are two different, related concepts:

- **Degree** is the number of edges touching a vertex.
- **Weighted Degree** is the sum of the weights of the edges touching a vertex.

For an unweighted graph, every edge has weight 1, so the two values are always equal, and only degree is shown. For a weighted graph, they can differ, so clicking a vertex shows both — for example, "Vertex 1 — Degree: 3, Weighted Degree: 4."

### Reading the Graph

The Graph Visualization card shows:

- **The graph itself** — vertices as numbered circles connected by straight edges, in a fixed, deterministic layout chosen from the graph's recognized shape (see below). The layout is the same every time for the same graph, so it never appears to "jump around." For a weighted graph, each edge's weight is drawn directly on it.
- **Vertices** and **Edges** counts.
- **Graph Type**, shown only when the shape can be identified with certainty, using standard mathematical notation with subscripts: **Complete Graph Kₙ**, **Path Graph Pₙ**, **Cycle Graph Cₙ**, **Star Graph K₁,ₙ₋₁**, **Disconnected Graph (X components)**, **Empty Graph (no edges)**, or **Single Vertex** — prefixed with "Weighted" for a weighted graph (for example, **Weighted Path Graph Pₙ**). If the graph doesn't match one of these recognizable shapes, an unweighted graph shows only the vertex/edge counts with no invented label, while a weighted graph still shows **Weighted Graph**, since the edge weights themselves remain useful information even without a recognized shape.
- A small **Weighted Graph** badge appears next to the Graph Visualization heading whenever the current graph is weighted, so it's clear at a glance without reading every edge.
- Disconnected graphs are drawn exactly as disconnected — the calculator never adds edges to force components together.

### Layout by Graph Shape

Rather than placing every graph on one generic circle, the calculator recognizes several common shapes directly from the reconstructed edges and draws each one the way it's normally drawn on paper. The matrix always determines the edges; the shape recognition below only ever changes *where* vertices are drawn, never *which* edges exist.

| Recognized shape | Layout used |
|---|---|
| **Path graph** (Pₙ) | A straight horizontal line, with vertices placed left-to-right in the order they actually connect along the path — not necessarily their numeric order. |
| **Cycle graph**, 4 vertices (C₄) | An axis-aligned square, one vertex per corner, in the order they connect around the cycle. |
| **Cycle graph**, other sizes (C₃, C₅, C₆, …) | A regular polygon (triangle, pentagon, hexagon, and so on), vertices placed around it in the order they connect around the cycle. |
| **Star graph** (K₁,ₙ₋₁) | The hub (the vertex connected to every other vertex) is placed dead center, with every leaf vertex arranged evenly around it. |
| **Complete graph**, 4 vertices (K₄) | A triangle with a central vertex: three vertices form the outer triangle and the fourth sits at the centroid, connected to all three — the standard way K₄ is drawn. |
| **Complete graph**, other sizes (K₃, K₅, K₆, …) | A symmetric regular polygon, one vertex per corner. |
| **Disconnected graph** | Each connected component is placed in its own clearly separated region of the drawing (arranged in a grid), and every component is independently laid out using the same shape rules above — so a disconnected graph made of two paths, for example, draws two separate lines rather than mixing every vertex onto one shared shape. |
| Anything else | A symmetric circular arrangement in vertex-number order, as a reliable general-purpose fallback. |

Path and cycle detection follows the graph's actual connections (not just vertex numbering), so the drawing stays a clean, uncrossed line or ring even if the underlying matrix numbers the vertices out of order along the path or cycle.

The table above describes the **unweighted** layout, where every edge is drawn the same length. A weighted graph keeps the same shape recognition (a weighted path is still drawn as a line, a weighted star still has a centered hub, and so on) but adapts each shape so edge length reflects weight, without changing which vertices are connected:

| Recognized shape (weighted) | How weight changes the layout |
|---|---|
| **Path graph** | Still a straight line, but each segment's length along that line is scaled by its edge's weight, so the shortest-weight segment is the shortest stretch of the line and the longest-weight segment is the longest. |
| **Cycle graph** | Still a ring of the same radius, but vertices are spaced unevenly around it — a lower-weight edge pulls its two vertices closer together along the ring, a higher-weight edge pushes them farther apart — so the ring stays a clean, uncrossed loop at any size while its arc lengths reflect weight. |
| **Star graph** | The hub still sits dead center, but each leaf sits at a distance from the hub proportional to that leaf's edge weight, instead of all leaves sitting on one fixed circle. |
| **Complete graph** (small) | Vertices are placed with a lightweight, bounded simulation — edges act like springs pulled toward a length proportional to their weight, while vertices gently repel each other so nothing overlaps — so a weighted K₄, for example, no longer has to be a perfect triangle-plus-center. |
| **Complete or unrecognized graph** (large) | Kept at the plain equal-length layout. Once a graph is large and dense enough that individual edges are already hard to trace (see [Large and Dense Graphs](#large-and-dense-graphs)), stretching edges by weight on top of that would only add visual noise, not clarity — the weight labels and hover tooltips remain the reliable way to read exact weights at that size. |

### Graph Interaction

- **Hover** (or keyboard-focus with Tab) a vertex to highlight it and every edge connected to it.
- **Hover** an edge, or its weight label, to highlight that edge on its own — its weight label is emphasized at the same time, and a tooltip shows the edge and its weight (for example, "Edge 1-2, Weight: 2").
- **Click** a vertex (or press Enter/Space while it's focused) to show its degree in a line below the graph — for example "Vertex 1 — degree 3" for an unweighted graph, or "Vertex 1 — Degree: 3, Weighted Degree: 4." for a weighted one.
- These interactions add convenience; none of them are required to read the graph — the vertex/edge counts, graph type, and the drawing itself (including edge weights) are always visible without hovering or clicking anything.

### Download Graph

Click **Download Graph** to save the current graph as a standalone `laplacian-graph.svg` file, containing only the graph drawing (vertices, labels, and edges) — not the rest of the page. For a weighted graph, the edge weight labels are included in the downloaded file exactly as shown on screen. SVG is a scalable, resolution-independent image format that any modern browser or image editor can open.

### Copy Graph Information

Click **Copy Graph Info** to copy a plain-text summary of the graph to your clipboard. For an unweighted graph:

```
Vertices: 4
Edges: 6
Graph Type: Complete Graph K₄

Edges:
1-2
1-3
1-4
2-3
2-4
3-4
```

For a weighted graph, each edge line also includes its weight:

```
Vertices: 4
Edges: 4
Graph Type: Weighted Graph

Edges:
1-2: weight 2
1-3: weight 1
1-4: weight 1
2-3: weight 3
```

### Large and Dense Graphs

The graph feature supports the calculator's full matrix size range, 1×1 through 100×100. As graphs get larger and denser:

- Vertex numbers are shown as permanent labels for smaller graphs, and on hover/focus only for graphs with more than roughly 45 vertices, to keep the drawing from becoming cluttered.
- Edge weight labels follow the same threshold: for a weighted graph with more than roughly 45 vertices, weights are available by hovering an edge (via the tooltip) rather than as permanent on-graph labels, again to avoid clutter.
- Edges are drawn with reduced opacity once a graph has many edges, so a very dense graph (for example, a complete graph with dozens of vertices) reads as a legible density pattern instead of a single solid, unreadable shape.
- A short note appears under large or dense graphs explaining what changed — for example, that individual edges may be hard to distinguish visually, or that edge weights are shown on hover instead of as permanent labels. The vertex and edge counts and the graph type (when identifiable) always remain accurate and readable regardless of density.
- The eigenvalue/eigenvector calculation itself is unaffected by graph size or density — the two features operate independently.

---

## Interface Guide

| Control | Purpose |
|---|---|
| Matrix Size (n) | Sets the dimension of the square matrix (1–100) |
| Create Matrix | Builds the matrix input grid for the chosen size |
| Calculate Eigenvalues & Eigenvectors | Computes and displays the eigenvalues and eigenvectors |
| Load Example | Loads a built-in 2×2 example matrix |
| Clear Matrix | Empties all cell values in the current matrix |
| Copy Results | Copies the matrix and calculated results to the clipboard |
| Copy Graph Info | Copies the vertex/edge summary of the reconstructed graph to the clipboard (including edge weights, for a weighted graph) |
| Download Graph | Downloads the current graph as an SVG image file (including edge weight labels, for a weighted graph) |

---

## Matrix Size and Large Matrices

The calculator supports square matrices from **1×1 up to 100×100**. Because a matrix has n × n entries, larger sizes mean significantly more cells to fill in:

- 10×10 = 100 entries
- 50×50 = 2,500 entries
- 100×100 = 10,000 entries

When you create a matrix of 50×50 or larger, the calculator displays a notice:

> "Large matrix detected. Calculations may take longer depending on your device."

Large matrices are still fully supported, but filling in thousands of cells by hand is naturally slower, and the calculation itself takes noticeably longer as the size grows. A 100×100 matrix typically calculates in a few seconds on a modern device — it is not instantaneous, and performance depends on your device's hardware.

---

## Technology Stack

- **HTML5** — page structure and layout
- **CSS3** — styling, responsive layout, and light/dark theming
- **JavaScript (vanilla, no frameworks)** — all interactivity, calculations, and graph rendering
- **Inline SVG** — the graph visualization, generated and drawn directly in the browser

No external math, linear algebra, or graph/charting library is used. The eigenvalue and eigenvector computations (Hessenberg reduction, the shifted QR algorithm with complex-number support, and null-space calculation) and the Laplacian validation, graph reconstruction, and SVG drawing are all implemented directly in the project's own JavaScript. There is no backend, server, or database — the entire calculator runs client-side in the browser.

---

## How It Works

1. The user selects a matrix size and clicks Create Matrix (or presses Enter).
2. The website generates an empty grid of input cells for that size.
3. The user enters a numeric value into every cell.
4. When Calculate is clicked, JavaScript reads and validates the matrix values.
5. The matrix is reduced to a simpler equivalent form, then processed with a numerical algorithm to find its eigenvalues.
6. For each eigenvalue found, the calculator determines its multiplicity and computes its associated eigenvector(s).
7. The eigenvalues, multiplicities, and eigenvectors are formatted and displayed on the page.
8. Independently of that calculation, the same matrix is checked against the Laplacian rules described above. If it qualifies, its vertices and edges (and, for a weighted Laplacian, each edge's weight) are reconstructed and drawn as an SVG graph; if not, an informational message explaining why is shown in its place.

All of this happens locally in the browser — no matrix data is ever transmitted anywhere.

---

## Project Structure

```text
Eigenvalue & Eigenvector Calculator/
├── .gitattributes
├── index.html
├── script.js
├── style.css
└── README.md
```

- **index.html** — the page structure and layout for the calculator
- **style.css** — all visual styling, responsive behavior, and theming
- **script.js** — the matrix interface logic and the eigenvalue/eigenvector calculation engine
- **.gitattributes** — line-ending normalization settings for the repository

---

## Running It Locally

This is a fully static website with no build step and no dependencies to install. To run it on your own machine:

1. Download or clone this repository.
2. Open the `Eigenvalue & Eigenvector Calculator` folder.
3. Open `index.html` directly in a web browser.

Alternatively, you can serve the folder with any simple local web server, for example:

```bash
# From inside the project folder
python -m http.server 8000
```

Then visit `http://localhost:8000` in your browser.

---

## Deployment

This project is a static site (HTML, CSS, and JavaScript only), so it can be deployed to any static hosting provider. It is currently deployed on **Netlify** at https://eigenvalueproj.netlify.app/.

To deploy your own copy on Netlify:

1. Push this repository to GitHub.
2. Create a new site on Netlify and connect it to the repository.
3. Since there is no build step, leave the build command empty and set the publish directory to the project folder (the one containing `index.html`).
4. Deploy — Netlify will serve the site directly.

The same approach works on any other static host (GitHub Pages, Vercel, Cloudflare Pages, etc.), since the project has no server-side requirements.

---

## Limitations & Notes

- This calculator uses floating-point numerical methods, not exact symbolic algebra. Results are typically accurate to several decimal places but may show very small rounding artifacts (for example, an extremely small value close to zero).
- Eigenvalues that are very close together (numerically clustered) can be harder to separate precisely, as with any numerical eigenvalue solver.
- Very large matrices (approaching 100×100) take noticeably longer to calculate and require entering a large number of values by hand.
- All matrix entries must be real numbers; the matrix itself must be square (the interface only allows square matrices, from 1×1 to 100×100).
- The calculator runs entirely in the browser using standard JavaScript, so performance depends on the device and browser being used.
- Graph visualization supports **simple** graphs, both unweighted and weighted: undirected, with no self-loops and no multiple edges between the same pair of vertices, matching the symmetric, zero-row-sum Laplacian rules described above. Directed graphs and multigraphs are not recognized as valid Laplacians and will not produce a graph.
- Edge weights must be positive numbers (a negative off-diagonal entry's absolute value). The calculator never invents or infers a weighted graph from a matrix that fails the Laplacian rules — if validation fails, no graph is shown, and the specific reason is explained instead.
- Graph type labeling (Complete, Path, Cycle, Star, Disconnected, Empty, Single Vertex) only covers those specific, unambiguous shapes. An unweighted graph that doesn't match one of them is still drawn correctly, just without a type label; a weighted graph that doesn't match one still shows the generic **Weighted Graph** label, since the weights themselves remain informative even without a recognized shape. The vertex and edge counts are always shown either way.
- The graph layout is a fixed circular arrangement; vertices cannot be manually dragged or rearranged. This keeps the layout deterministic and reliable at every size, including very large graphs, at the cost of custom positioning.
- For graphs with a large number of edges, individual edges can be visually difficult to trace by eye even though every edge is drawn correctly — the vertex/edge counts and graph type remain accurate regardless.
