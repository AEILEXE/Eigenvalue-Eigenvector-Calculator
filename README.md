# Eigenvalue & Eigenvector Calculator

A browser-based calculator that computes the eigenvalues and eigenvectors of any square matrix. Enter a matrix size, fill in the values, and get the eigenvalues, their multiplicities, and their corresponding eigenvectors — all calculated instantly in your browser, with no data ever sent to a server. When the matrix you enter is a graph Laplacian, the calculator also reconstructs and displays the underlying graph automatically.

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
- Automatically visualize the graph behind a Laplacian matrix, when the matrix you entered is one

In linear algebra, an eigenvalue and eigenvector pair describes a direction that a matrix only stretches or shrinks, without changing — a property used throughout engineering, physics, computer graphics, and data science (for example, in vibration analysis, stability analysis, and principal component analysis). This calculator handles that computation for you, including matrices whose eigenvalues turn out to be complex numbers.

A graph Laplacian is a special square matrix (`L = D - A`, degree matrix minus adjacency matrix) built from a graph's vertices and edges; its eigenvalues reveal structural properties of the graph (for example, the number of zero eigenvalues equals the number of connected components). This calculator recognizes valid graph Laplacians automatically and draws the graph they represent alongside the usual eigenvalue results.

---

## Features

- Square matrix input for any size from 1×1 to 100×100
- Eigenvalue calculation, including complex (non-real) eigenvalues
- Eigenvector calculation for every eigenvalue
- Algebraic multiplicity (how many times an eigenvalue repeats as a root)
- Geometric multiplicity (how many independent eigenvectors an eigenvalue actually has)
- A "Defective" indicator when an eigenvalue's algebraic and geometric multiplicities differ
- Automatic Laplacian detection and graph reconstruction, with a clean, interactive SVG visualization
- Automatic graph type labeling for recognizable shapes (complete, path, cycle, star, disconnected, or empty graphs)
- Hover and keyboard-focus highlighting of graph vertices and their edges, and click-to-reveal vertex degree
- Downloadable graph image (SVG) and one-click copying of the graph's vertex/edge summary
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

If the matrix you calculated is a valid simple-graph Laplacian, a **Graph Visualization** card automatically appears below the eigenvalue results, showing the graph the matrix represents. If it isn't, the card instead shows a short message explaining that graph visualization isn't available for that matrix — the eigenvalue and eigenvector results above are unaffected either way. See [Laplacian Matrix and Graph Visualization](#laplacian-matrix-and-graph-visualization) below for full details.

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

Whenever you calculate a matrix that fits this pattern, the calculator reconstructs the graph directly from the matrix entries and displays it automatically — you don't need to do anything extra.

### How the Calculator Recognizes a Laplacian

A calculated matrix is treated as a valid simple-graph Laplacian only when **all** of the following hold:

- It is square (guaranteed by the matrix size you chose).
- Every entry is a finite real number.
- It is symmetric (`L[i][j] = L[j][i]` for every pair).
- Every row sums to zero (within a small numerical tolerance).
- Every diagonal entry is a non-negative whole number.
- Every off-diagonal entry is either `0` or `-1`.
- Each diagonal entry equals the count of `-1` entries in that same row (its vertex's degree).

If any of these checks fail, the matrix is not treated as a Laplacian. The eigenvalue and eigenvector results are calculated and shown exactly as normal, and the Graph Visualization card simply displays:

> "Graph visualization is unavailable because this matrix is not a valid simple graph Laplacian."

This is an informational message, not an error — the rest of the calculator continues to work normally.

### How the Graph Is Reconstructed

The graph is built directly from the matrix, never from the eigenvalues:

- Each row/column index becomes a vertex, numbered starting from **1** (row/column 1 in the matrix is vertex 1, and so on).
- For every off-diagonal entry `L[i][j] = -1` (with i < j), one edge is drawn between vertex i and vertex j. Because the matrix is symmetric, the mirrored entry `L[j][i] = -1` is not read again, so each edge is created exactly once.

For example, the Laplacian

```
 3 -1 -1 -1
-1  3 -1 -1
-1 -1  3 -1
-1 -1 -1  3
```

reconstructs into 4 vertices, each connected to every other vertex — 6 edges in total (the complete graph K₄).

### Reading the Graph

The Graph Visualization card shows:

- **The graph itself** — vertices as numbered circles connected by straight edges, in a fixed, deterministic layout chosen from the graph's recognized shape (see below). The layout is the same every time for the same graph, so it never appears to "jump around."
- **Vertices** and **Edges** counts.
- **Graph Type**, shown only when the shape can be identified with certainty, using standard mathematical notation with subscripts: **Complete Graph Kₙ**, **Path Graph Pₙ**, **Cycle Graph Cₙ**, **Star Graph K₁,ₙ₋₁**, **Disconnected Graph (X components)**, **Empty Graph (no edges)**, or **Single Vertex**. If the graph doesn't match one of these recognizable shapes, only the vertex/edge counts are shown, with no invented label.
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

### Graph Interaction

- **Hover** (or keyboard-focus with Tab) a vertex to highlight it and every edge connected to it.
- **Hover** an edge to highlight that edge on its own.
- **Click** a vertex (or press Enter/Space while it's focused) to show its degree in a line below the graph, for example "Vertex 1 — degree 3."
- These interactions add convenience; none of them are required to read the graph — the vertex/edge counts, graph type, and the drawing itself are always visible without hovering or clicking anything.

### Download Graph

Click **Download Graph** to save the current graph as a standalone `laplacian-graph.svg` file, containing only the graph drawing (vertices, labels, and edges) — not the rest of the page. SVG is a scalable, resolution-independent image format that any modern browser or image editor can open.

### Copy Graph Information

Click **Copy Graph Info** to copy a plain-text summary of the graph to your clipboard, for example:

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

### Large and Dense Graphs

The graph feature supports the calculator's full matrix size range, 1×1 through 100×100. As graphs get larger and denser:

- Vertex numbers are shown as permanent labels for smaller graphs, and on hover/focus only for graphs with more than roughly 45 vertices, to keep the drawing from becoming cluttered.
- Edges are drawn with reduced opacity once a graph has many edges, so a very dense graph (for example, a complete graph with dozens of vertices) reads as a legible density pattern instead of a single solid, unreadable shape.
- A short note appears under large or dense graphs noting that individual edges may be hard to distinguish visually. The vertex and edge counts and the graph type (when identifiable) always remain accurate and readable regardless of density.
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
| Copy Graph Info | Copies the vertex/edge summary of the reconstructed graph to the clipboard |
| Download Graph | Downloads the current graph as an SVG image file |

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
8. Independently of that calculation, the same matrix is checked against the Laplacian rules described above. If it qualifies, its vertices and edges are reconstructed and drawn as an SVG graph; if not, an informational message is shown in its place.

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
- Graph visualization only supports **simple** graphs: unweighted, undirected, with no self-loops and no multiple edges between the same pair of vertices, matching the strict `0`/`-1` off-diagonal rule described above. Weighted Laplacians, directed graphs, and multigraphs are not recognized as valid Laplacians and will not produce a graph.
- Graph type labeling (Complete, Path, Cycle, Star, Disconnected, Empty, Single Vertex) only covers those specific, unambiguous shapes. A valid graph that doesn't match one of them is still drawn correctly, just without a type label — the vertex and edge counts are always shown either way.
- The graph layout is a fixed circular arrangement; vertices cannot be manually dragged or rearranged. This keeps the layout deterministic and reliable at every size, including very large graphs, at the cost of custom positioning.
- For graphs with a large number of edges, individual edges can be visually difficult to trace by eye even though every edge is drawn correctly — the vertex/edge counts and graph type remain accurate regardless.
