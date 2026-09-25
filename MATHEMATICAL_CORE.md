# Mathematical Core

## 1. Purpose and Scope

This document is the mathematical and algorithmic reference for the Eigenvalue & Eigenvector Calculator. It covers exactly two things: (1) the numerical linear-algebra engine that computes eigenvalues, eigenvectors, and their multiplicities, and (2) the graph mathematics derived from a matrix recognized as a Laplacian (edge/weight extraction, degree, weighted degree, connected components, and shape classification).

**Out of scope:** HTML markup, CSS, color/typography/theming, page layout, button and card styling, animations, decorative SVG attributes, graph *layout/positioning* (pixel coordinates), DOM event wiring, and any other presentation-layer detail. Those live in `index.html`, `style.css`, and the UI portions of `script.js`; they are not discussed here. For a complete description of the system as a whole — UI, workflow, features — see `README.md`, which is the system-level documentation and defers to this document for mathematical detail.

Everything quoted in this document is copied from the current `script.js`, not paraphrased, simplified, or reconstructed from memory. Every numeric example in [Section 14](#14-important-mathematical-examples) was produced by extracting the exact functions below out of the current `script.js` and executing them in Node.js against the stated input matrix — not computed by hand or estimated.

---

## 2. Mathematical Pipeline

Two independent pipelines both start from the same validated `n×n` array of real numbers. Neither reads the other's output.

**Eigenvalue / eigenvector pipeline:**

```
Matrix
  → Matrix Processing        (realToComplexMatrix, toHessenberg)
  → Eigenvalue Calculation   (eigenvaluesQR: complexGivens + hessenbergQRStep, shifted QR with deflation)
  → Eigenvalue Grouping      (groupEigenvalues: clusters numerically-equal roots)
  → Algebraic Multiplicity   (= cluster size, from groupEigenvalues)
  → Nullspace                (nullspace(matrix, n, eigenvalue), per distinct eigenvalue)
  → Eigenvectors             (= the nullspace basis)
  → Geometric Multiplicity   (= basis.length, i.e. vectors.length)
  → Defective Detection      (vectors.length < multiplicity)
```

**Laplacian / graph-data pipeline:**

```
Matrix
  → Laplacian Validation      (laplacianAnalysis: symmetry, sign, row-sum checks)
  → Edge Extraction           (same function: L[i][j] < 0  ->  edge (i, j), weight = |L[i][j]|)
  → Degree Calculation        (computeDegrees)
  → Weighted Degree           (computeWeightedDegrees, only when weighted)
  → Connected Components      (countComponents)
  → Graph Classification      (classifyGraph: Kn / Pn / Cn / star / disconnected / empty / single)
```

---

## 3. Complex Arithmetic

The eigenvalue solver runs entirely in complex arithmetic so that a complex-conjugate eigenvalue pair is produced directly by the same code path as a real eigenvalue, with no separate case. Every function below is a primitive that the rest of the pipeline is built from.

```js
function cx(re, im) { return { re: re, im: im || 0 }; }
function cAdd(a, b) { return cx(a.re + b.re, a.im + b.im); }
function cSub(a, b) { return cx(a.re - b.re, a.im - b.im); }
function cMul(a, b) { return cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
function cDiv(a, b) {
  var d = b.re * b.re + b.im * b.im;
  if (d === 0) d = 1e-300;
  return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
function cAbs(a) { return Math.hypot(a.re, a.im); }
function cNeg(a) { return cx(-a.re, -a.im); }
function cConj(a) { return cx(a.re, -a.im); }
function cScaleR(a, s) { return cx(a.re * s, a.im * s); }
function cSqrt(z) {
  var r = cAbs(z);
  if (r === 0) return cx(0, 0);
  var re = Math.sqrt((r + z.re) / 2);
  var im = Math.sqrt(Math.max(0, r - z.re) / 2);
  if (z.im < 0) im = -im;
  return cx(re, im);
}
```
*(`script.js:27-47`)*

| Function | Role |
|---|---|
| `cx` | Constructs a `{ re, im }` complex number (imaginary part defaults to 0). |
| `cAdd` / `cSub` / `cMul` / `cDiv` | The four arithmetic operations, applied component-wise per the standard complex-number formulas. |
| `cAbs` | Magnitude `\|z\|`, via `Math.hypot` for numerical stability. |
| `cNeg` / `cConj` | Negation and complex conjugate. |
| `cScaleR` | Multiplies a complex number by a real scalar. |
| `cSqrt` | Principal complex square root (needed for the quadratic-formula shift computation in Section 5, since the discriminant can be negative or complex). |

---

## 4. Matrix Processing

The real input matrix is converted to complex form, then reduced to upper-Hessenberg form before the QR algorithm runs. Hessenberg reduction (via real Householder reflections) is a **similarity transform**, so it changes the matrix's shape but not its eigenvalues — it exists purely to make each QR iteration in Section 5 cheap (`O(m²)` per step instead of `O(m³)`).

```js
/* ----------------------------- Complex matrices ----------------------------- */
function realToComplexMatrix(A, n) {
  return A.map(function (row) { return row.map(function (v) { return cx(v, 0); }); });
}

/* ----------------------------- Hessenberg reduction ----------------------------- */
function toHessenberg(realA, n) {
  var A = realA.map(function (row) { return row.slice(); });
  for (var k = 0; k < n - 2; k++) {
    var L = n - k - 1;
    var x = [];
    for (var i = k + 1; i < n; i++) x.push(A[i][k]);
    var normx = Math.sqrt(x.reduce(function (s, v) { return s + v * v; }, 0));
    if (normx < 1e-300) continue;
    var alpha = x[0] >= 0 ? -normx : normx;
    var v = x.slice();
    v[0] -= alpha;
    var vnorm2 = v.reduce(function (s, vv) { return s + vv * vv; }, 0);
    if (vnorm2 < 1e-300) continue;
    var beta = 2 / vnorm2;

    for (var j = 0; j < n; j++) {
      var s = 0;
      for (var l = 0; l < L; l++) s += v[l] * A[k + 1 + l][j];
      s *= beta;
      for (l = 0; l < L; l++) A[k + 1 + l][j] -= v[l] * s;
    }
    for (var r = 0; r < n; r++) {
      var s2 = 0;
      for (l = 0; l < L; l++) s2 += A[r][k + 1 + l] * v[l];
      s2 *= beta;
      for (l = 0; l < L; l++) A[r][k + 1 + l] -= s2 * v[l];
    }
  }
  return A;
}
```
*(`script.js:50-84`)*

For each column `k`, `toHessenberg` builds a Householder reflector `v` that zeroes every entry below the subdiagonal in that column, then applies it from both sides (`A ← (I − βvvᵀ)A(I − βvvᵀ)`) so the result stays similar to the original matrix. After `n − 2` columns, `A` is upper-Hessenberg: zero everywhere below the first subdiagonal.

---

## 5. QR Eigenvalue Algorithm

This is the actual implementation, not a generic textbook description of QR iteration.

`complexGivens` computes one Givens rotation; `hessenbergQRStep` applies `m − 1` of them to perform one **implicit-shift QR step** on an `m × m` upper-Hessenberg matrix. Because a QR step on a Hessenberg matrix produces another Hessenberg matrix, this is mathematically the same shifted-QR iteration a dense implementation would perform — only the per-step cost changes (`O(m²)` via Givens rotations instead of `O(m³)` via a dense Householder QR decomposition + matrix multiply). This is what keeps a 100×100 matrix practical to run in a browser tab.

```js
/* ----------------------------- Hessenberg-preserving QR step -----------------------------
   One implicit-shift QR step on an m x m upper-Hessenberg matrix, done with
   m-1 Givens rotations (O(m^2)) instead of a dense Householder QR
   decomposition + matrix multiply (O(m^3)). A QR step on a Hessenberg
   matrix produces another Hessenberg matrix, so this is mathematically
   the same shifted-QR iteration as a dense implementation -- only the
   cost of computing each step changes, not what is computed. This keeps
   large matrices (up to 100x100) practical: a dense-QR-per-step approach
   took ~16.5s end-to-end for a 100x100 matrix; this reduces that to ~3s. */
function complexGivens(a, b) {
  if (cAbs(b) < 1e-300) return { c: 1, s: cx(0, 0) };
  if (cAbs(a) < 1e-300) return { c: 0, s: cx(1, 0) };
  var absa = cAbs(a);
  var norm = Math.sqrt(absa * absa + cAbs(b) * cAbs(b));
  var c = absa / norm;
  var phaseA = cx(a.re / absa, a.im / absa);
  var s = cScaleR(cMul(phaseA, cConj(b)), 1 / norm);
  return { c: c, s: s };
}

function hessenbergQRStep(M, m, shift) {
  var i, j, k, row, col;
  for (i = 0; i < m; i++) M[i][i] = cSub(M[i][i], shift);

  var rotations = [];
  for (k = 0; k < m - 1; k++) {
    var g = complexGivens(M[k][k], M[k + 1][k]);
    rotations.push(g);
    for (col = k; col < m; col++) {
      var top = M[k][col], bot = M[k + 1][col];
      M[k][col] = cAdd(cMul(cx(g.c, 0), top), cMul(g.s, bot));
      M[k + 1][col] = cAdd(cMul(cNeg(cConj(g.s)), top), cMul(cx(g.c, 0), bot));
    }
    M[k + 1][k] = cx(0, 0);
  }

  for (k = 0; k < m - 1; k++) {
    var gg = rotations[k];
    for (row = 0; row < m; row++) {
      var left = M[row][k], right = M[row][k + 1];
      M[row][k] = cAdd(cMul(left, cx(gg.c, 0)), cMul(right, cConj(gg.s)));
      M[row][k + 1] = cAdd(cMul(left, cNeg(gg.s)), cMul(right, cx(gg.c, 0)));
    }
  }

  for (i = 0; i < m; i++) M[i][i] = cAdd(M[i][i], shift);
}

/* ----------------------------- Shifted QR eigenvalue solver ----------------------------- */
function eigenvaluesQR(realA, n) {
  if (n === 1) return [cx(realA[0][0], 0)];
  var hess = toHessenberg(realA, n);
  var A = realToComplexMatrix(hess, n);
  var eigenvalues = [];
  var m = n;
  var maxTotalIter = 100 * n + 500;
  var iterCount = 0;

  while (m > 1 && iterCount < maxTotalIter) {
    var stagnant = 0;
    var prevSub = Infinity;
    while (m > 1) {
      iterCount++;
      if (iterCount > maxTotalIter) break;
      var sub = cAbs(A[m - 1][m - 2]);
      var scale = cAbs(A[m - 2][m - 2]) + cAbs(A[m - 1][m - 1]);
      if (scale < 1e-300) scale = 1;
      if (sub < 1e-13 * scale) break;

      var a = A[m - 2][m - 2], b = A[m - 2][m - 1], c = A[m - 1][m - 2], d = A[m - 1][m - 1];
      var tr = cAdd(a, d);
      var det = cSub(cMul(a, d), cMul(b, c));
      var disc = cSqrt(cSub(cMul(tr, tr), cScaleR(det, 4)));
      var e1 = cScaleR(cAdd(tr, disc), 0.5);
      var e2 = cScaleR(cSub(tr, disc), 0.5);
      var shift = cAbs(cSub(e1, d)) < cAbs(cSub(e2, d)) ? e1 : e2;

      if (Math.abs(sub - prevSub) < 1e-15 * scale) {
        stagnant++;
        if (stagnant > 10) {
          shift = cAdd(shift, cx(scale * 0.1, scale * 0.13));
          stagnant = 0;
        }
      }
      prevSub = sub;

      var subM = [];
      for (var i = 0; i < m; i++) subM.push(A[i].slice(0, m));
      hessenbergQRStep(subM, m, shift);
      for (i = 0; i < m; i++) for (var j = 0; j < m; j++) A[i][j] = subM[i][j];
    }
    eigenvalues.push(A[m - 1][m - 1]);
    m--;
  }
  if (m === 1) eigenvalues.push(A[0][0]);
  return eigenvalues;
}
```
*(`script.js:95-182`)*

**Shift selection:** each iteration solves the 2×2 eigenvalue problem for the trailing `2×2` block (`a, b, c, d`) directly via the quadratic formula (`tr`, `det`, `disc = √(tr² − 4·det)`, giving candidate shifts `e1`, `e2`), then picks whichever of the two is closer to `d` — this is the Wilkinson-style shift strategy, chosen because it converges quickly and handles complex eigenvalue pairs without a separate real/complex branch.

**Stagnation escape:** if the subdiagonal entry being driven toward zero (`sub`) stops changing between iterations (within `1e-15 * scale`) for more than 10 consecutive iterations, the shift is perturbed by a small fixed complex offset (`scale * 0.1 + i·scale * 0.13`) to break the stall.

**Deflation:** once the subdiagonal entry `A[m-1][m-2]` is smaller than `1e-13 * scale`, it is treated as converged; the trailing eigenvalue `A[m-1][m-1]` is recorded and the working matrix size `m` is decremented by one (`m--`), so subsequent iterations work on a strictly smaller trailing submatrix.

**Iteration limit:** `maxTotalIter = 100 * n + 500` bounds the total number of QR steps across the whole run, as a safety limit against non-convergence; the algorithm returns whatever eigenvalues have been deflated so far if that limit is hit.

**Base case:** for `n === 1`, the single matrix entry is returned directly as a real-valued complex number — no Hessenberg reduction or QR iteration is needed.

---

## 6. Eigenvalue Grouping

`eigenvaluesQR` returns exactly `n` roots. For a genuinely repeated eigenvalue, floating-point arithmetic returns several numerically-close-but-not-bit-identical values rather than one exact duplicate. `groupEigenvalues` clusters roots that are within a scale-relative tolerance of each other and averages each cluster into a single reported value — **the size of a cluster is the algebraic multiplicity**.

```js
/* ----------------------------- Grouping (algebraic multiplicity) ----------------------------- */
function internalClean(x, scale) {
  var tol = 1e-6 * Math.max(1, scale || 1);
  if (Math.abs(x) < tol) return 0;
  var r = Math.round(x);
  if (Math.abs(x - r) < tol) return r;
  return Math.round(x * 1e6) / 1e6;
}

function groupEigenvalues(roots) {
  var scale = 1;
  roots.forEach(function (r) { scale = Math.max(scale, Math.abs(r.re), Math.abs(r.im)); });
  var clusterTol = Math.max(1e-4, 1e-5 * scale);
  var used = new Array(roots.length).fill(false);
  var groups = [];
  for (var i = 0; i < roots.length; i++) {
    if (used[i]) continue;
    var members = [roots[i]];
    used[i] = true;
    for (var j = i + 1; j < roots.length; j++) {
      if (used[j]) continue;
      var dist = Math.hypot(roots[i].re - roots[j].re, roots[i].im - roots[j].im);
      if (dist < clusterTol) { members.push(roots[j]); used[j] = true; }
    }
    var avgRe = members.reduce(function (s, m) { return s + m.re; }, 0) / members.length;
    var avgIm = members.reduce(function (s, m) { return s + m.im; }, 0) / members.length;
    groups.push({
      value: cx(internalClean(avgRe, scale), internalClean(avgIm, scale)),
      rawValue: cx(avgRe, avgIm),
      multiplicity: members.length
    });
  }
  groups.sort(function (a, b) {
    if (a.value.re !== b.value.re) return a.value.re - b.value.re;
    return a.value.im - b.value.im;
  });
  return groups;
}
```
*(`script.js:185-221`)*

- **Clustering tolerance** (`clusterTol`) is `max(1e-4, 1e-5 * scale)`, where `scale` is the largest real or imaginary magnitude across all roots — so the tolerance scales with the matrix's own numbers rather than using one fixed absolute value that would be too loose for tiny eigenvalues or too tight for huge ones.
- Clustering is a simple greedy pairwise pass: pick the first unused root, absorb every remaining unused root within `clusterTol` of it (by Euclidean distance in the complex plane), then move to the next unused root.
- Each group's reported `value` is the cluster's average, passed through `internalClean` to snap near-zero values to exactly `0` and near-integer values to the exact integer (both within a `1e-6 * scale` tolerance) — this is what makes an eigenvalue like `3.0000002` display and compare as `3`. The unrounded average is kept separately as `rawValue`, which is what actually gets passed into `nullspace` (Section 7) for numerical accuracy.
- **`members.length` is the algebraic multiplicity**, stored as `multiplicity` on each group.
- Groups are sorted ascending, first by real part, then by imaginary part.

---

## 7. Nullspace and Eigenvectors

For a given eigenvalue λ, its eigenvectors are exactly the null space of `(A − λI)` — i.e. every nonzero vector `v` satisfying:

```
(A − λI) v = 0
```

`nullspace` finds a basis for that null space using complex Gaussian elimination with partial pivoting (row-reduction to a reduced-row-echelon-like form), then reads off one basis vector per free (non-pivot) column.

```js
/* ----------------------------- Eigenvector null space ----------------------------- */
function nullspace(A, n, lambda) {
  var B = [];
  for (var i = 0; i < n; i++) {
    B.push([]);
    for (var j = 0; j < n; j++) {
      var v = cx(A[i][j], 0);
      if (i === j) v = cSub(v, lambda);
      B[i].push(v);
    }
  }
  var maxEntry = 0;
  for (i = 0; i < n; i++) for (j = 0; j < n; j++) maxEntry = Math.max(maxEntry, cAbs(B[i][j]));
  var eps = 1e-6 * Math.max(1, maxEntry);

  var pivotRow = 0;
  var pivotCols = [];
  var pivotRowForCol = {};
  for (var col = 0; col < n && pivotRow < n; col++) {
    var maxVal = eps;
    var maxR = -1;
    for (var r = pivotRow; r < n; r++) {
      var m = cAbs(B[r][col]);
      if (m > maxVal) { maxVal = m; maxR = r; }
    }
    if (maxR === -1) continue;
    var tmp = B[pivotRow]; B[pivotRow] = B[maxR]; B[maxR] = tmp;

    var pivotVal = B[pivotRow][col];
    for (var c2 = 0; c2 < n; c2++) B[pivotRow][c2] = cDiv(B[pivotRow][c2], pivotVal);

    for (var r2 = 0; r2 < n; r2++) {
      if (r2 === pivotRow) continue;
      var factor = B[r2][col];
      if (cAbs(factor) < 1e-300) continue;
      for (var c3 = 0; c3 < n; c3++) {
        B[r2][c3] = cSub(B[r2][c3], cMul(factor, B[pivotRow][c3]));
      }
    }
    pivotCols.push(col);
    pivotRowForCol[col] = pivotRow;
    pivotRow++;
  }

  var freeCols = [];
  for (col = 0; col < n; col++) if (pivotCols.indexOf(col) === -1) freeCols.push(col);

  var basis = [];
  freeCols.forEach(function (f) {
    var v = new Array(n).fill(null).map(function () { return cx(0, 0); });
    v[f] = cx(1, 0);
    pivotCols.forEach(function (pc) {
      var row = pivotRowForCol[pc];
      v[pc] = cNeg(B[row][f]);
    });
    basis.push(v);
  });
  return basis;
}
```
*(`script.js:224-281`)*

- `B` is built as `(A − λI)` in complex form: `λ` is subtracted only from the diagonal.
- `eps = 1e-6 * max(1, maxEntry)` is the pivoting tolerance — scaled to the matrix's own largest entry magnitude, again to stay meaningful across very different matrix scales.
- The elimination loop is standard partial-pivoting Gaussian elimination adapted to complex numbers: for each column, find the row (at or below the current pivot row) with the largest-magnitude entry in that column; if it clears `eps`, swap it into the pivot position, scale that row so the pivot becomes `1`, then eliminate that column's entry from every *other* row (this is why the result is a reduced form, not just triangular — it makes reading off the basis in the next step a direct lookup rather than a back-substitution).
- Columns that never receive a pivot are the **free columns** — each one corresponds to one dimension of the null space, i.e. one eigenvector.
- The basis is constructed directly: for each free column `f`, set that vector's `f`-th entry to `1` and every pivot-column entry to `−B[row][f]` (the negated reduced coefficient) — this is the standard "set free variable to 1, solve for the pivot variables" construction, done once per free column.

---

## 8. Geometric Multiplicity

There is no dedicated "geometric multiplicity" function. Geometric multiplicity is, by definition, **the dimension of the eigenspace** for a given eigenvalue — and since `nullspace` already returns a basis for exactly that eigenspace, geometric multiplicity is simply **the length of the array `nullspace` returns** (`freeCols.length`, i.e. `basis.length`).

`computeEigen` is the function that assembles the full pipeline from Sections 5–7: it computes the roots, groups them into algebraic multiplicities, and then calls `nullspace` once per distinct eigenvalue to attach that eigenvalue's eigenvector basis.

```js
/* ----------------------------- Public entry point ----------------------------- */
function computeEigen(matrix, n) {
  if (n === 1) {
    var val = matrix[0][0];
    return [{ value: cx(internalClean(val, Math.abs(val) || 1), 0), multiplicity: 1, vectors: [[cx(1, 0)]] }];
  }
  var roots = eigenvaluesQR(matrix, n);
  var groups = groupEigenvalues(roots);
  groups.forEach(function (g) {
    g.vectors = nullspace(matrix, n, g.rawValue);
  });
  return groups;
}
```
*(`script.js:284-295`)*

For each returned group `g`:

- `g.multiplicity` — the **algebraic** multiplicity (from Section 6: how many roots clustered into this eigenvalue).
- `g.vectors.length` — the **geometric** multiplicity (from Section 7: the dimension of the null space of `A − λI`).

---

## 9. Defective Matrices

An eigenvalue is classified as **defective** wherever its eigenspace has fewer dimensions than its multiplicity as a root — i.e. exactly when:

```
geometric multiplicity < algebraic multiplicity
```
which, in terms of the values computed above, is the literal comparison:
```js
g.vectors.length < g.multiplicity
```
This is exactly the condition the rest of the application uses to decide whether to show the "Defective" indicator for a given eigenvalue result. There is no separate defectiveness-detection function — it is this one comparison, evaluated per eigenvalue group after Sections 6–7 have both run. A defective eigenvalue means the matrix does not have a full set of `n` linearly independent eigenvectors overall, i.e. `A` is not diagonalizable.

---

## 10. Laplacian Mathematics

For a simple, undirected graph, the graph Laplacian is the matrix:

```
L = D − A
```

where `D` is the diagonal degree matrix and `A` is the adjacency matrix. `laplacianAnalysis` is the single function that validates whether an arbitrary user-entered matrix is a valid Laplacian — weighted or unweighted; an unweighted Laplacian is simply the special case where every edge weight equals exactly `1`. It checks, **in this exact order**, stopping and returning a specific reason at the first failure:

1. Every entry is a finite real number.
2. The matrix is symmetric: `L[i][j] = L[j][i]` for every pair.
3. Every off-diagonal entry is zero or negative (a positive off-diagonal entry is never valid for a Laplacian).
4. Every diagonal entry is non-negative.
5. Every row sums to zero, within a numerical tolerance.

No condition beyond these five is tested — in particular, there is no separate check that the matrix is square (it always is, by construction, since the UI only ever builds `n×n` grids), and no check on the size or magnitude of individual weights beyond what row-sum-zero and the sign rule already imply.

```js
var LAPLACIAN_TOL = 1e-6;

/* A matrix L is a valid (weighted or unweighted) simple-graph Laplacian
   (L = D - A) when it is symmetric, every off-diagonal entry is zero or
   negative (a negative entry's magnitude is the edge weight), every
   diagonal entry is non-negative, and every row sums to zero. An
   unweighted Laplacian is simply the special case where every edge
   weight equals 1 -- so this single analysis handles both, and the
   original K4-style 0/-1 matrices validate exactly as before. Returns
   either { ok: true, edges: [[i, j, weight], ...], weighted: bool } or
   { ok: false, reason: <human-readable explanation> } -- never silently
   "fixes" or reinterprets an invalid matrix. */
function laplacianAnalysis(matrix, n) {
  if (!n || n < 1) return { ok: false, reason: 'the matrix is empty.' };
  var i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (!isFinite(matrix[i][j])) return { ok: false, reason: 'the matrix contains a non-finite value.' };
    }
  }
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      if (Math.abs(matrix[i][j] - matrix[j][i]) > LAPLACIAN_TOL) {
        return {
          ok: false,
          reason: 'the matrix is not symmetric (L[' + (i + 1) + '][' + (j + 1) + '] = ' + formatReal(matrix[i][j]) +
            ' but L[' + (j + 1) + '][' + (i + 1) + '] = ' + formatReal(matrix[j][i]) + '), so it cannot represent an undirected graph.'
        };
      }
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i === j) continue;
      if (matrix[i][j] > LAPLACIAN_TOL) {
        return {
          ok: false,
          reason: 'entry L[' + (i + 1) + '][' + (j + 1) + '] = ' + formatReal(matrix[i][j]) +
            ' is positive; off-diagonal entries of a Laplacian must be zero or negative (a negative value represents an edge weight).'
        };
      }
    }
    if (matrix[i][i] < -LAPLACIAN_TOL) {
      return {
        ok: false,
        reason: 'diagonal entry L[' + (i + 1) + '][' + (i + 1) + '] is negative; Laplacian diagonal entries must be non-negative.'
      };
    }
  }
  for (i = 0; i < n; i++) {
    var rowSum = 0;
    for (j = 0; j < n; j++) rowSum += matrix[i][j];
    if (Math.abs(rowSum) > LAPLACIAN_TOL * Math.max(1, n)) {
      return {
        ok: false,
        reason: 'row ' + (i + 1) + ' sums to ' + formatReal(rowSum) + ', not 0; every row of a Laplacian matrix must sum to zero.'
      };
    }
  }

  var edges = [];
  var weighted = false;
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      var val = matrix[i][j];
      if (val < -LAPLACIAN_TOL) {
        var w = -val;
        edges.push([i + 1, j + 1, w]);
        if (Math.abs(w - 1) > LAPLACIAN_TOL) weighted = true;
      }
    }
  }
  return { ok: true, edges: edges, weighted: weighted };
}
```
*(`LAPLACIAN_TOL` at `script.js:695`; function body at `script.js:708-769`.)*

All tolerance comparisons use the single constant `LAPLACIAN_TOL = 1e-6`, except the row-sum check, which scales it by `n` (`LAPLACIAN_TOL * Math.max(1, n)`) so the acceptable floating-point drift grows with how many terms are being summed. If validation fails at any step, the function returns `{ ok: false, reason }` with a specific human-readable explanation of exactly which rule failed and why — the matrix is never silently corrected or reinterpreted. If every check passes, it returns `{ ok: true, edges, weighted }`.

---

## 11. Graph Data Mathematics

This section covers the data derived from an already-validated Laplacian: edges, weights, degree, weighted degree, and connected components. None of it involves pixels, coordinates, or SVG.

**Edge and weight extraction** happens in the same pass at the end of `laplacianAnalysis` (Section 10). This is the mathematically important distinction the rest of this document (and `README.md`) relies on:

> **MATHEMATICAL EDGE WEIGHT** is `w = |L[i][j]|`, derived directly from the Laplacian, exactly as computed below.
> It is a completely different thing from **VISUAL EDGE LENGTH** — the number of on-screen pixels an edge is drawn with. Edge weight is graph data; edge length is a rendering decision made later, entirely outside this document's scope (in `script.js`'s layout code, not shown here — see [Section 12](#12-weighted-graph-mathematics) for exactly where the boundary is).

```js
var edges = [];
var weighted = false;
for (i = 0; i < n; i++) {
  for (j = i + 1; j < n; j++) {
    var val = matrix[i][j];
    if (val < -LAPLACIAN_TOL) {
      var w = -val;
      edges.push([i + 1, j + 1, w]);
      if (Math.abs(w - 1) > LAPLACIAN_TOL) weighted = true;
    }
  }
}
return { ok: true, edges: edges, weighted: weighted };
```
*(the closing slice of `laplacianAnalysis`, `script.js:756-768` — already shown in full in Section 10; repeated here as the specific excerpt this section is about)*

Only the upper triangle (`j > i`) is scanned, since the matrix is already known symmetric at this point — so each edge is discovered exactly once, never twice. Each edge is stored as a 3-element array `[vertex_i, vertex_j, weight]`, with vertex numbers 1-indexed (row/column `0` in the matrix is vertex `1`). The graph is flagged `weighted` the moment any single edge's weight differs from `1` by more than the tolerance — one non-unit-weight edge is enough to make the whole graph "weighted."

**Degree and weighted degree**, computed directly from that edge list:

```js
function computeDegrees(n, edges) {
  var deg = new Array(n).fill(0);
  edges.forEach(function (e) {
    deg[e[0] - 1]++;
    deg[e[1] - 1]++;
  });
  return deg;
}

/* Weighted degree: sum of incident edge weights (distinct from plain
   degree, which is the count of incident edges). */
function computeWeightedDegrees(n, edges) {
  var wd = new Array(n).fill(0);
  edges.forEach(function (e) {
    wd[e[0] - 1] += e[2];
    wd[e[1] - 1] += e[2];
  });
  return wd;
}
```
*(`script.js:771-789`)*

`computeDegrees` counts incident edges per vertex (each edge increments both of its endpoints' counters). `computeWeightedDegrees` sums incident edge *weights* per vertex instead of counting edges — for an unweighted graph (every weight `1`) the two functions necessarily produce identical results, which is why the UI only bothers computing/showing weighted degree when `weighted` is true. Note that `L[i][i]` (the diagonal) is *never* read directly as the weighted degree anywhere in the code — the diagonal is only ever used to *validate* the row-sum-zero rule in Section 10; weighted degree is always recomputed independently from the edge list, as shown here.

**Connected components**, a union-find (disjoint-set-union) over the vertex set, used by graph classification to detect a disconnected graph:

```js
function countComponents(n, edges) {
  var parent = [];
  for (var i = 0; i <= n; i++) parent[i] = i;
  function find(x) { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
  edges.forEach(function (e) {
    var ra = find(e[0]), rb = find(e[1]);
    if (ra !== rb) parent[ra] = rb;
  });
  var roots = {};
  for (var v = 1; v <= n; v++) roots[find(v)] = true;
  return Object.keys(roots).length;
}
```
*(`script.js:791-802`)*

Each vertex starts as its own parent (its own component); for every edge, the two endpoints' components are unioned (with path compression in `find` for efficiency). The number of distinct roots left at the end is the number of connected components.

---

## 12. Weighted Graph Mathematics

For a Laplacian matrix `L`, the weight of the edge between vertex `i` and vertex `j` (when that edge exists) is derived directly from the negative off-diagonal entry:

```
w_ij = −L_ij            (for i ≠ j, whenever L_ij < 0)
```

`L_ij = 0` means no edge exists between `i` and `j` — it is not interpreted as "an edge with weight 0." This is exactly the rule implemented in the `edges.push([i + 1, j + 1, w])` line quoted in full in Section 11: `val` is `L_ij`, and `w = -val` is the weight, computed only inside the `if (val < -LAPLACIAN_TOL)` branch (i.e. only when an edge actually exists).

An **unweighted** Laplacian is the special case of this same rule where every existing edge's weight happens to equal `1` — `laplacianAnalysis` does not run separate code paths for the weighted and unweighted cases; it is one analysis, and the `weighted` boolean it returns is just a summary of whether any extracted weight differs from `1`.

**Weighted degree** (already shown in full in Section 11, `computeWeightedDegrees`) is the sum, per vertex, of every incident edge's weight `w_ij` — mathematically:

```
weighted_degree(v) = Σ w_vj   over every edge (v, j)
```

This is distinct from plain **degree**, `Σ 1` over the same edges (i.e. just the edge count) — the two are computed by separate functions (`computeDegrees` and `computeWeightedDegrees`) because they answer different questions and are numerically identical only in the unweighted case.

**Where this document's scope ends:** the weight `w_ij` computed above is graph *data*. What happens to that number afterward — normalizing it against the graph's minimum/maximum weight and mapping it to a pixel length so smaller-weight edges draw visually shorter — is a rendering/layout decision made in a separate part of `script.js` (functions such as `desiredEdgeLength`, `computeWeightedLayoutPoints`, and the force-directed/analytic layout routines for paths, cycles, and stars). None of that code is reproduced in this document: it consumes the mathematical data above but does not compute anything mathematically new about the graph itself, and per this document's scope (Section 1), presentation/layout code is intentionally excluded. `README.md` documents that visual behavior at the system level.

---

## 13. Graph Classification

`classifyGraph` assigns a structural label to a graph purely from its vertex count, edge count, and per-vertex degree sequence — it never looks at weight, so a weighted and unweighted graph with the same underlying connections classify identically (the "Weighted" prefix seen in the UI is a separate string concatenation done outside this function, not part of the classification logic itself).

```js
var SUBSCRIPT_DIGITS = { '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉' };
function toSubscript(num) {
  return String(num).split('').map(function (ch) { return SUBSCRIPT_DIGITS[ch] || ch; }).join('');
}

/* Only classifies a graph shape when the check is unambiguous; otherwise
   returns null so the UI falls back to a plain vertex/edge count. Uses
   standard mathematical notation (Kn, Pn, Cn, K1,n-1) with Unicode
   subscripts for the vertex count. */
function classifyGraph(n, edges, degrees) {
  var m = edges.length;
  if (n === 1) return 'Single Vertex';
  if (m === 0) return 'Empty Graph (no edges)';

  var components = countComponents(n, edges);
  if (components > 1) return 'Disconnected Graph (' + components + ' components)';

  var maxEdges = (n * (n - 1)) / 2;
  if (m === maxEdges) return 'Complete Graph K' + toSubscript(n);

  var deg1Count = 0, deg2Count = 0, otherDeg = false, maxDeg = 0;
  degrees.forEach(function (d) {
    if (d === 1) deg1Count++;
    else if (d === 2) deg2Count++;
    else otherDeg = true;
    if (d > maxDeg) maxDeg = d;
  });

  if (!otherDeg && m === n - 1 && deg1Count === 2 && deg2Count === n - 2) {
    return 'Path Graph P' + toSubscript(n);
  }
  if (!otherDeg && m === n && deg2Count === n && n >= 3) {
    return 'Cycle Graph C' + toSubscript(n);
  }
  if (m === n - 1 && maxDeg === n - 1 && deg1Count === n - 1) {
    return 'Star Graph K' + toSubscript(1) + ',' + toSubscript(n - 1);
  }
  return null;
}
```
*(`script.js:804-842`)*

The checks are evaluated in this order, each one final (no re-checking once a label is chosen):

1. `n === 1` → **Single Vertex**.
2. Zero edges → **Empty Graph (no edges)**.
3. More than one connected component (via `countComponents`, Section 11) → **Disconnected Graph (X components)** — checked *before* any shape test, so a disconnected graph is never mistaken for a connected shape.
4. Edge count equals `n(n−1)/2`, the maximum possible for `n` vertices → **Complete Graph Kₙ**.
5. Degree sequence is exactly two vertices of degree 1 and the rest degree 2, with `m = n − 1` edges → **Path Graph Pₙ**.
6. Every vertex has degree exactly 2 and `m = n` (and `n ≥ 3`) → **Cycle Graph Cₙ**.
7. One vertex has degree `n − 1` (connected to everyone) and every other vertex has degree exactly 1, with `m = n − 1` edges → **Star Graph K₁,ₙ₋₁**.
8. None of the above match → returns `null` (no label; the UI shows only the vertex/edge counts, or, for a weighted graph, a generic "Weighted Graph" label — that fallback wording is UI logic, not part of this function).

---

## 14. Important Mathematical Examples

Every result below was produced by extracting the exact functions from Sections 3–13 out of the current `script.js` and executing them directly in Node.js against the stated matrix — not computed by hand, not estimated, and not carried over from an earlier/different version of the code.

### Eigenvalue / eigenvector examples

| Case | Matrix | Verified eigenvalues (λ, algebraic mult., geometric mult.) |
|---|---|---|
| 1×1 | `[7]` | `7` (1, 1) |
| Distinct (2×2, "Load Example") | `[[2,1],[0,3]]` | `2` (1, 1); `3` (1, 1) |
| Repeated, non-defective (2×2) | `[[2,0],[0,2]]` | `2` (2, 2) |
| **Defective** (2×2, Jordan block) | `[[2,1],[0,2]]` | `2` (2, **1** — DEFECTIVE) |
| Complex (2×2, rotation) | `[[0,-1],[1,0]]` | `−i` (1, 1); `i` (1, 1) |
| Identity (3×3) | `[[1,0,0],[0,1,0],[0,0,1]]` | `1` (3, 3) |
| Zero matrix (3×3) | `[[0,0,0],[0,0,0],[0,0,0]]` | `0` (3, 3) |
| Diagonal (3×3) | `[[5,0,0],[0,-2,0],[0,0,7]]` | `−2` (1, 1); `5` (1, 1); `7` (1, 1) |

### Graph Laplacian examples (unweighted)

| Graph | Matrix | Verified eigenvalues | Verified classification |
|---|---|---|---|
| K₃ | `[[2,-1,-1],[-1,2,-1],[-1,-1,2]]` | `0` (1,1); `3` (2,2) | vertices=3, edges=3, `Complete Graph K₃` |
| P₄ | `[[1,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,1]]` | `0`; `0.585786`; `2`; `3.414214` (all mult. 1) | vertices=4, edges=3, `Path Graph P₄` |
| K₁,₃ | `[[3,-1,-1,-1],[-1,1,0,0],[-1,0,1,0],[-1,0,0,1]]` | `0` (1,1); `1` (2,2); `4` (1,1) | vertices=4, edges=3, `Star Graph K₁,₃` |
| C₄ | `[[2,-1,0,-1],[-1,2,-1,0],[0,-1,2,-1],[-1,0,-1,2]]` | `0` (1,1); `2` (2,2); `4` (1,1) | vertices=4, edges=4, `Cycle Graph C₄` |
| K₄ | `[[3,-1,-1,-1],[-1,3,-1,-1],[-1,-1,3,-1],[-1,-1,-1,3]]` | `0` (1,1); `4` (3,3) | vertices=4, edges=6, `Complete Graph K₄` |
| Disconnected (two separate edges, 1–2 and 3–4) | `[[1,-1,0,0],[-1,1,0,0],[0,0,1,-1],[0,0,-1,1]]` | `0` (2,2); `2` (2,2) | vertices=4, edges=2, `Disconnected Graph (2 components)` |

### Weighted graph Laplacian examples

**Weighted K₄** — this is the canonical verified example for this project. Matrix:

```
 7  -2  -1  -4
-2   8  -3  -3
-1  -3   5  -1
-4  -3  -1   8
```

Extracted edges (weight = `|L_ij|`): 1–2 = 2, 1–3 = 1, 1–4 = 4, 2–3 = 3, 2–4 = 3, 3–4 = 1. Verified classification: vertices=4, edges=6, `weighted=true`, `Complete Graph K₄` (rendered by the UI as "Weighted Complete Graph K₄").

**Verified eigenvalues:**

```
0
5.763932
10.236068
12
```

each with algebraic multiplicity 1 and geometric multiplicity 1 (none repeated, none defective).

> **Important:** `0, 4, 8, 16` is **not** a correct eigenvalue set for this matrix — `det(A − 4I)`, `det(A − 8I)`, and `det(A − 16I)` are all nonzero. The verified set above (`0, 5.763932, 10.236068, 12`) is the one produced by the current `eigenvaluesQR`/`groupEigenvalues` implementation and independently confirmed via direct determinant evaluation at each candidate value.

**Weighted P₄.** Matrix:

```
 2  -2   0   0
-2   5  -3   0
 0  -3   4  -1
 0   0  -1   1
```

Edges: 1–2 = 2, 2–3 = 3, 3–4 = 1. Verified classification: vertices=4, edges=3, `weighted=true`, `Path Graph P₄`. Verified eigenvalues: `0`, `1`, `3`, `8` (each multiplicity 1).

**Weighted C₄.** Matrix:

```
 4  -1   0  -3
-1   5  -4   0
 0  -4   6  -2
-3   0  -2   5
```

Edges: 1–2 = 1, 1–4 = 3, 2–3 = 4, 3–4 = 2. Verified classification: vertices=4, edges=4, `weighted=true`, `Cycle Graph C₄`. Verified eigenvalues: `0`, `2.855148`, `6.718356`, `10.426496` (each multiplicity 1).

**Weighted star** (hub = vertex 1; leaves 2, 3, 4 with weights 2, 4, 1 respectively). Matrix:

```
 7  -2  -4  -1
-2   2   0   0
-4   0   4   0
-1   0   0   1
```

Verified classification: vertices=4, edges=3, `weighted=true`, `Star Graph K₁,₃`. Verified eigenvalues: `0`, `1.202521`, `2.612831`, `10.184648` (each multiplicity 1).

### Invalid-matrix examples (verified rejection reasons)

| Matrix | `laplacianAnalysis` result |
|---|---|
| `[[3,-2,-1],[-2,5,0],[-1,0,1]]` | rejected: *"row 2 sums to 3, not 0; every row of a Laplacian matrix must sum to zero."* |
| `[[2,-2,0],[-1,3,-2],[0,-2,2]]` | rejected: *"the matrix is not symmetric (L[1][2] = -2 but L[2][1] = -1), so it cannot represent an undirected graph."* |
| `[[2,1],[0,3]]` (the "Load Example" matrix) | rejected: *"the matrix is not symmetric (L[1][2] = 1 but L[2][1] = 0), so it cannot represent an undirected graph."* — confirming that an ordinary (non-Laplacian) matrix still runs through `computeEigen` normally; it simply produces no graph. |

---

# Complete Mathematical Core Source Code

```javascript
function cx(re, im) { return { re: re, im: im || 0 }; }
function cAdd(a, b) { return cx(a.re + b.re, a.im + b.im); }
function cSub(a, b) { return cx(a.re - b.re, a.im - b.im); }
function cMul(a, b) { return cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
function cDiv(a, b) {
  var d = b.re * b.re + b.im * b.im;
  if (d === 0) d = 1e-300;
  return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
function cAbs(a) { return Math.hypot(a.re, a.im); }
function cNeg(a) { return cx(-a.re, -a.im); }
function cConj(a) { return cx(a.re, -a.im); }
function cScaleR(a, s) { return cx(a.re * s, a.im * s); }
function cSqrt(z) {
  var r = cAbs(z);
  if (r === 0) return cx(0, 0);
  var re = Math.sqrt((r + z.re) / 2);
  var im = Math.sqrt(Math.max(0, r - z.re) / 2);
  if (z.im < 0) im = -im;
  return cx(re, im);
}

/* ----------------------------- Complex matrices ----------------------------- */
function realToComplexMatrix(A, n) {
  return A.map(function (row) { return row.map(function (v) { return cx(v, 0); }); });
}

/* ----------------------------- Hessenberg reduction ----------------------------- */
function toHessenberg(realA, n) {
  var A = realA.map(function (row) { return row.slice(); });
  for (var k = 0; k < n - 2; k++) {
    var L = n - k - 1;
    var x = [];
    for (var i = k + 1; i < n; i++) x.push(A[i][k]);
    var normx = Math.sqrt(x.reduce(function (s, v) { return s + v * v; }, 0));
    if (normx < 1e-300) continue;
    var alpha = x[0] >= 0 ? -normx : normx;
    var v = x.slice();
    v[0] -= alpha;
    var vnorm2 = v.reduce(function (s, vv) { return s + vv * vv; }, 0);
    if (vnorm2 < 1e-300) continue;
    var beta = 2 / vnorm2;

    for (var j = 0; j < n; j++) {
      var s = 0;
      for (var l = 0; l < L; l++) s += v[l] * A[k + 1 + l][j];
      s *= beta;
      for (l = 0; l < L; l++) A[k + 1 + l][j] -= v[l] * s;
    }
    for (var r = 0; r < n; r++) {
      var s2 = 0;
      for (l = 0; l < L; l++) s2 += A[r][k + 1 + l] * v[l];
      s2 *= beta;
      for (l = 0; l < L; l++) A[r][k + 1 + l] -= s2 * v[l];
    }
  }
  return A;
}

/* ----------------------------- Hessenberg-preserving QR step -----------------------------
   One implicit-shift QR step on an m x m upper-Hessenberg matrix, done with
   m-1 Givens rotations (O(m^2)) instead of a dense Householder QR
   decomposition + matrix multiply (O(m^3)). A QR step on a Hessenberg
   matrix produces another Hessenberg matrix, so this is mathematically
   the same shifted-QR iteration as a dense implementation -- only the
   cost of computing each step changes, not what is computed. This keeps
   large matrices (up to 100x100) practical: a dense-QR-per-step approach
   took ~16.5s end-to-end for a 100x100 matrix; this reduces that to ~3s. */
function complexGivens(a, b) {
  if (cAbs(b) < 1e-300) return { c: 1, s: cx(0, 0) };
  if (cAbs(a) < 1e-300) return { c: 0, s: cx(1, 0) };
  var absa = cAbs(a);
  var norm = Math.sqrt(absa * absa + cAbs(b) * cAbs(b));
  var c = absa / norm;
  var phaseA = cx(a.re / absa, a.im / absa);
  var s = cScaleR(cMul(phaseA, cConj(b)), 1 / norm);
  return { c: c, s: s };
}

function hessenbergQRStep(M, m, shift) {
  var i, j, k, row, col;
  for (i = 0; i < m; i++) M[i][i] = cSub(M[i][i], shift);

  var rotations = [];
  for (k = 0; k < m - 1; k++) {
    var g = complexGivens(M[k][k], M[k + 1][k]);
    rotations.push(g);
    for (col = k; col < m; col++) {
      var top = M[k][col], bot = M[k + 1][col];
      M[k][col] = cAdd(cMul(cx(g.c, 0), top), cMul(g.s, bot));
      M[k + 1][col] = cAdd(cMul(cNeg(cConj(g.s)), top), cMul(cx(g.c, 0), bot));
    }
    M[k + 1][k] = cx(0, 0);
  }

  for (k = 0; k < m - 1; k++) {
    var gg = rotations[k];
    for (row = 0; row < m; row++) {
      var left = M[row][k], right = M[row][k + 1];
      M[row][k] = cAdd(cMul(left, cx(gg.c, 0)), cMul(right, cConj(gg.s)));
      M[row][k + 1] = cAdd(cMul(left, cNeg(gg.s)), cMul(right, cx(gg.c, 0)));
    }
  }

  for (i = 0; i < m; i++) M[i][i] = cAdd(M[i][i], shift);
}

/* ----------------------------- Shifted QR eigenvalue solver ----------------------------- */
function eigenvaluesQR(realA, n) {
  if (n === 1) return [cx(realA[0][0], 0)];
  var hess = toHessenberg(realA, n);
  var A = realToComplexMatrix(hess, n);
  var eigenvalues = [];
  var m = n;
  var maxTotalIter = 100 * n + 500;
  var iterCount = 0;

  while (m > 1 && iterCount < maxTotalIter) {
    var stagnant = 0;
    var prevSub = Infinity;
    while (m > 1) {
      iterCount++;
      if (iterCount > maxTotalIter) break;
      var sub = cAbs(A[m - 1][m - 2]);
      var scale = cAbs(A[m - 2][m - 2]) + cAbs(A[m - 1][m - 1]);
      if (scale < 1e-300) scale = 1;
      if (sub < 1e-13 * scale) break;

      var a = A[m - 2][m - 2], b = A[m - 2][m - 1], c = A[m - 1][m - 2], d = A[m - 1][m - 1];
      var tr = cAdd(a, d);
      var det = cSub(cMul(a, d), cMul(b, c));
      var disc = cSqrt(cSub(cMul(tr, tr), cScaleR(det, 4)));
      var e1 = cScaleR(cAdd(tr, disc), 0.5);
      var e2 = cScaleR(cSub(tr, disc), 0.5);
      var shift = cAbs(cSub(e1, d)) < cAbs(cSub(e2, d)) ? e1 : e2;

      if (Math.abs(sub - prevSub) < 1e-15 * scale) {
        stagnant++;
        if (stagnant > 10) {
          shift = cAdd(shift, cx(scale * 0.1, scale * 0.13));
          stagnant = 0;
        }
      }
      prevSub = sub;

      var subM = [];
      for (var i = 0; i < m; i++) subM.push(A[i].slice(0, m));
      hessenbergQRStep(subM, m, shift);
      for (i = 0; i < m; i++) for (var j = 0; j < m; j++) A[i][j] = subM[i][j];
    }
    eigenvalues.push(A[m - 1][m - 1]);
    m--;
  }
  if (m === 1) eigenvalues.push(A[0][0]);
  return eigenvalues;
}

/* ----------------------------- Grouping (algebraic multiplicity) ----------------------------- */
function internalClean(x, scale) {
  var tol = 1e-6 * Math.max(1, scale || 1);
  if (Math.abs(x) < tol) return 0;
  var r = Math.round(x);
  if (Math.abs(x - r) < tol) return r;
  return Math.round(x * 1e6) / 1e6;
}

function groupEigenvalues(roots) {
  var scale = 1;
  roots.forEach(function (r) { scale = Math.max(scale, Math.abs(r.re), Math.abs(r.im)); });
  var clusterTol = Math.max(1e-4, 1e-5 * scale);
  var used = new Array(roots.length).fill(false);
  var groups = [];
  for (var i = 0; i < roots.length; i++) {
    if (used[i]) continue;
    var members = [roots[i]];
    used[i] = true;
    for (var j = i + 1; j < roots.length; j++) {
      if (used[j]) continue;
      var dist = Math.hypot(roots[i].re - roots[j].re, roots[i].im - roots[j].im);
      if (dist < clusterTol) { members.push(roots[j]); used[j] = true; }
    }
    var avgRe = members.reduce(function (s, m) { return s + m.re; }, 0) / members.length;
    var avgIm = members.reduce(function (s, m) { return s + m.im; }, 0) / members.length;
    groups.push({
      value: cx(internalClean(avgRe, scale), internalClean(avgIm, scale)),
      rawValue: cx(avgRe, avgIm),
      multiplicity: members.length
    });
  }
  groups.sort(function (a, b) {
    if (a.value.re !== b.value.re) return a.value.re - b.value.re;
    return a.value.im - b.value.im;
  });
  return groups;
}

/* ----------------------------- Eigenvector null space ----------------------------- */
function nullspace(A, n, lambda) {
  var B = [];
  for (var i = 0; i < n; i++) {
    B.push([]);
    for (var j = 0; j < n; j++) {
      var v = cx(A[i][j], 0);
      if (i === j) v = cSub(v, lambda);
      B[i].push(v);
    }
  }
  var maxEntry = 0;
  for (i = 0; i < n; i++) for (j = 0; j < n; j++) maxEntry = Math.max(maxEntry, cAbs(B[i][j]));
  var eps = 1e-6 * Math.max(1, maxEntry);

  var pivotRow = 0;
  var pivotCols = [];
  var pivotRowForCol = {};
  for (var col = 0; col < n && pivotRow < n; col++) {
    var maxVal = eps;
    var maxR = -1;
    for (var r = pivotRow; r < n; r++) {
      var m = cAbs(B[r][col]);
      if (m > maxVal) { maxVal = m; maxR = r; }
    }
    if (maxR === -1) continue;
    var tmp = B[pivotRow]; B[pivotRow] = B[maxR]; B[maxR] = tmp;

    var pivotVal = B[pivotRow][col];
    for (var c2 = 0; c2 < n; c2++) B[pivotRow][c2] = cDiv(B[pivotRow][c2], pivotVal);

    for (var r2 = 0; r2 < n; r2++) {
      if (r2 === pivotRow) continue;
      var factor = B[r2][col];
      if (cAbs(factor) < 1e-300) continue;
      for (var c3 = 0; c3 < n; c3++) {
        B[r2][c3] = cSub(B[r2][c3], cMul(factor, B[pivotRow][c3]));
      }
    }
    pivotCols.push(col);
    pivotRowForCol[col] = pivotRow;
    pivotRow++;
  }

  var freeCols = [];
  for (col = 0; col < n; col++) if (pivotCols.indexOf(col) === -1) freeCols.push(col);

  var basis = [];
  freeCols.forEach(function (f) {
    var v = new Array(n).fill(null).map(function () { return cx(0, 0); });
    v[f] = cx(1, 0);
    pivotCols.forEach(function (pc) {
      var row = pivotRowForCol[pc];
      v[pc] = cNeg(B[row][f]);
    });
    basis.push(v);
  });
  return basis;
}

/* ----------------------------- Public entry point ----------------------------- */
function computeEigen(matrix, n) {
  if (n === 1) {
    var val = matrix[0][0];
    return [{ value: cx(internalClean(val, Math.abs(val) || 1), 0), multiplicity: 1, vectors: [[cx(1, 0)]] }];
  }
  var roots = eigenvaluesQR(matrix, n);
  var groups = groupEigenvalues(roots);
  groups.forEach(function (g) {
    g.vectors = nullspace(matrix, n, g.rawValue);
  });
  return groups;
}

function formatReal(x) {
  if (Math.abs(x) < 1e-9) x = 0;
  var rounded = Math.round(x * 1e6) / 1e6;
  if (Object.is(rounded, -0)) rounded = 0;
  var str = rounded.toFixed(6);
  str = str.replace(/(\.\d*?)0+$/, '$1').replace(/\.$/, '');
  if (str === '' || str === '-0') str = '0';
  return str;
}

var LAPLACIAN_TOL = 1e-6;

/* A matrix L is a valid (weighted or unweighted) simple-graph Laplacian
   (L = D - A) when it is symmetric, every off-diagonal entry is zero or
   negative (a negative entry's magnitude is the edge weight), every
   diagonal entry is non-negative, and every row sums to zero. An
   unweighted Laplacian is simply the special case where every edge
   weight equals 1 -- so this single analysis handles both, and the
   original K4-style 0/-1 matrices validate exactly as before. Returns
   either { ok: true, edges: [[i, j, weight], ...], weighted: bool } or
   { ok: false, reason: <human-readable explanation> } -- never silently
   "fixes" or reinterprets an invalid matrix. */
function laplacianAnalysis(matrix, n) {
  if (!n || n < 1) return { ok: false, reason: 'the matrix is empty.' };
  var i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (!isFinite(matrix[i][j])) return { ok: false, reason: 'the matrix contains a non-finite value.' };
    }
  }
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      if (Math.abs(matrix[i][j] - matrix[j][i]) > LAPLACIAN_TOL) {
        return {
          ok: false,
          reason: 'the matrix is not symmetric (L[' + (i + 1) + '][' + (j + 1) + '] = ' + formatReal(matrix[i][j]) +
            ' but L[' + (j + 1) + '][' + (i + 1) + '] = ' + formatReal(matrix[j][i]) + '), so it cannot represent an undirected graph.'
        };
      }
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i === j) continue;
      if (matrix[i][j] > LAPLACIAN_TOL) {
        return {
          ok: false,
          reason: 'entry L[' + (i + 1) + '][' + (j + 1) + '] = ' + formatReal(matrix[i][j]) +
            ' is positive; off-diagonal entries of a Laplacian must be zero or negative (a negative value represents an edge weight).'
        };
      }
    }
    if (matrix[i][i] < -LAPLACIAN_TOL) {
      return {
        ok: false,
        reason: 'diagonal entry L[' + (i + 1) + '][' + (i + 1) + '] is negative; Laplacian diagonal entries must be non-negative.'
      };
    }
  }
  for (i = 0; i < n; i++) {
    var rowSum = 0;
    for (j = 0; j < n; j++) rowSum += matrix[i][j];
    if (Math.abs(rowSum) > LAPLACIAN_TOL * Math.max(1, n)) {
      return {
        ok: false,
        reason: 'row ' + (i + 1) + ' sums to ' + formatReal(rowSum) + ', not 0; every row of a Laplacian matrix must sum to zero.'
      };
    }
  }

  var edges = [];
  var weighted = false;
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      var val = matrix[i][j];
      if (val < -LAPLACIAN_TOL) {
        var w = -val;
        edges.push([i + 1, j + 1, w]);
        if (Math.abs(w - 1) > LAPLACIAN_TOL) weighted = true;
      }
    }
  }
  return { ok: true, edges: edges, weighted: weighted };
}

function computeDegrees(n, edges) {
  var deg = new Array(n).fill(0);
  edges.forEach(function (e) {
    deg[e[0] - 1]++;
    deg[e[1] - 1]++;
  });
  return deg;
}

/* Weighted degree: sum of incident edge weights (distinct from plain
   degree, which is the count of incident edges). */
function computeWeightedDegrees(n, edges) {
  var wd = new Array(n).fill(0);
  edges.forEach(function (e) {
    wd[e[0] - 1] += e[2];
    wd[e[1] - 1] += e[2];
  });
  return wd;
}

function countComponents(n, edges) {
  var parent = [];
  for (var i = 0; i <= n; i++) parent[i] = i;
  function find(x) { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
  edges.forEach(function (e) {
    var ra = find(e[0]), rb = find(e[1]);
    if (ra !== rb) parent[ra] = rb;
  });
  var roots = {};
  for (var v = 1; v <= n; v++) roots[find(v)] = true;
  return Object.keys(roots).length;
}

var SUBSCRIPT_DIGITS = { '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉' };
function toSubscript(num) {
  return String(num).split('').map(function (ch) { return SUBSCRIPT_DIGITS[ch] || ch; }).join('');
}

/* Only classifies a graph shape when the check is unambiguous; otherwise
   returns null so the UI falls back to a plain vertex/edge count. Uses
   standard mathematical notation (Kn, Pn, Cn, K1,n-1) with Unicode
   subscripts for the vertex count. */
function classifyGraph(n, edges, degrees) {
  var m = edges.length;
  if (n === 1) return 'Single Vertex';
  if (m === 0) return 'Empty Graph (no edges)';

  var components = countComponents(n, edges);
  if (components > 1) return 'Disconnected Graph (' + components + ' components)';

  var maxEdges = (n * (n - 1)) / 2;
  if (m === maxEdges) return 'Complete Graph K' + toSubscript(n);

  var deg1Count = 0, deg2Count = 0, otherDeg = false, maxDeg = 0;
  degrees.forEach(function (d) {
    if (d === 1) deg1Count++;
    else if (d === 2) deg2Count++;
    else otherDeg = true;
    if (d > maxDeg) maxDeg = d;
  });

  if (!otherDeg && m === n - 1 && deg1Count === 2 && deg2Count === n - 2) {
    return 'Path Graph P' + toSubscript(n);
  }
  if (!otherDeg && m === n && deg2Count === n && n >= 3) {
    return 'Cycle Graph C' + toSubscript(n);
  }
  if (m === n - 1 && maxDeg === n - 1 && deg1Count === n - 1) {
    return 'Star Graph K' + toSubscript(1) + ',' + toSubscript(n - 1);
  }
  return null;
}
```
