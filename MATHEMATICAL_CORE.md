# Mathematical Core

## 1. Purpose

This document isolates the actual computational "brain" of the current website — the JavaScript functions in `script.js` that perform the linear algebra (eigenvalues, eigenvectors, multiplicities) and the graph/Laplacian data processing (validation, weighted/unweighted detection, edge and weight extraction, degree calculation).

Everything below is copied verbatim from the current `script.js`, as it exists after the weighted-graph revision. Nothing here is invented, simplified into pseudocode, or backported from an earlier version. Line numbers refer to the current `script.js`.

This file intentionally excludes HTML, CSS, SVG drawing, graph *layout/positioning* code, and DOM/UI event handling — see [Section 12](#12-what-is-not-part-of-the-mathematical-core) for exactly what was left out and why.

---

## 2. Main Mathematical Pipeline

There are two independent pipelines that both start from the same validated numeric matrix. Neither one depends on the other's output — the eigenvalue engine never reads graph data, and the graph analysis never reads eigenvalues.

**Eigenvalue / eigenvector pipeline:**

```
matrix (n×n array of real numbers)
  → computeEigen(matrix, n)
      → eigenvaluesQR(matrix, n)            // toHessenberg + shifted QR (Givens rotations)
      → groupEigenvalues(roots)             // clusters equal roots -> algebraic multiplicity
      → nullspace(matrix, n, eigenvalue)     // per group -> eigenvectors, whose count is geometric multiplicity
  → [{ value, multiplicity, vectors }, ...]
```

**Laplacian / graph data pipeline:**

```
matrix (the same n×n array)
  → laplacianAnalysis(matrix, n)
      → { ok: false, reason }                          // invalid Laplacian: no graph, eigen pipeline above is unaffected
      → { ok: true, edges: [[i, j, weight], ...], weighted }
  → computeDegrees(n, edges)               // plain degree
  → computeWeightedDegrees(n, edges)       // weighted degree (only when weighted)
  → classifyGraph(n, edges, degrees)       // Kn / Pn / Cn / star / disconnected / etc.
```

The actual orchestration that wires both pipelines to the matrix a user entered:

```js
function runCalculation() {
  showError('');
  var matrix = readMatrix();
  if (!matrix) return;

  calculateBtn.disabled = true;
  calculatingIndicator.hidden = false;

  window.setTimeout(function () {
    try {
      var groups = computeEigen(matrix, currentN);
      renderResults(matrix, currentN, groups);
      renderGraphSection(matrix, currentN);
    } catch (err) {
      showError('Something went wrong while calculating. Please check your matrix values and try again.');
    } finally {
      calculateBtn.disabled = false;
      calculatingIndicator.hidden = true;
    }
  }, 30);
}
```
*(`script.js:1681-1701`. `readMatrix()`, `renderResults()`, and `renderGraphSection()`'s DOM work are UI, not shown here — this excerpt exists only to show where `computeEigen` and the Laplacian pipeline are actually called from.)*

The first few lines of `renderGraphSection` show exactly how the Laplacian pipeline's output feeds into the graph data functions, before any DOM/SVG rendering happens:

```js
function renderGraphSection(matrix, n) {
  var analysis = laplacianAnalysis(matrix, n);
  if (!analysis.ok) {
    // ... DOM: shows analysis.reason in the "unavailable" message, no graph is built ...
    lastGraph = null;
    return;
  }

  // ... DOM: unhide the graph card ...

  var edges = analysis.edges;
  var weighted = analysis.weighted;
  var degrees = computeDegrees(n, edges);
  var weightedDegrees = weighted ? computeWeightedDegrees(n, edges) : null;
  var svg = buildGraphSVG(n, edges, degrees, weighted, weightedDegrees);
  // ... DOM/SVG rendering, classifyGraph() for the type label, etc. (not shown here) ...
}
```
*(`script.js:1555-1574`, trimmed — the omitted lines are `innerHTML`/`hidden`/`aria-label` DOM writes, not math.)*

---

## 3. Matrix Processing

Converts the real input matrix into the complex representation and upper-Hessenberg form the QR solver needs. Hessenberg reduction is a similarity transform (via real Householder reflections), so it does not change the matrix's eigenvalues — it only makes the shifted-QR iteration in Section 5 cheap to run.

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

---

## 4. Complex Number Operations

The QR eigenvalue solver runs entirely in complex arithmetic (so it can return complex eigenvalue pairs directly, without a separate special case). These are the primitive operations everything else in Sections 5–7 is built from.

```js
/* ----------------------------- Complex numbers ----------------------------- */
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

---

## 5. Eigenvalue Calculation

The shifted-QR algorithm with Wilkinson shifts, run on the Hessenberg matrix from Section 3. `complexGivens`/`hessenbergQRStep` perform one QR step via `m-1` Givens rotations instead of a full dense QR decomposition (mathematically identical result, much cheaper per step — this is what keeps a 100×100 matrix practical in-browser). `eigenvaluesQR` repeatedly applies this step with deflation until every eigenvalue (real or complex) has been extracted.

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

---

## 6. Eigenvalue Grouping and Algebraic Multiplicity

`eigenvaluesQR` returns `n` roots, which for a repeated eigenvalue arrive as several numerically-close-but-not-identical complex values (floating point). `groupEigenvalues` clusters roots within a relative tolerance and averages each cluster — the cluster size *is* the algebraic multiplicity.

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
*(`script.js:185-221`. `members.length` on line 213, stored as `multiplicity`, is the algebraic multiplicity.)*

---

## 7. Eigenvector Calculation

For a given eigenvalue λ, the eigenvectors are a basis for the null space of `(A − λI)`, found by complex Gaussian elimination (reduced row-echelon form) followed by reading off one basis vector per free column.

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

---

## 8. Geometric Multiplicity

There is no separate "geometric multiplicity" function — geometric multiplicity is simply **the number of basis vectors `nullspace()` returns** for a given eigenvalue (`freeCols.length`, i.e. `basis.length` in Section 7). `computeEigen` is the function that ties everything in Sections 5–7 together: it calls `eigenvaluesQR`, groups the roots into algebraic multiplicities, and then calls `nullspace` once per distinct eigenvalue to attach eigenvectors — the count of which the rest of the app reads as `g.vectors.length` to determine geometric multiplicity, and to flag an eigenvalue as **defective** when `g.vectors.length < g.multiplicity`.

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
*(`script.js:284-295`. `g.multiplicity` = algebraic multiplicity (from Section 6); `g.vectors.length` = geometric multiplicity (from Section 7).)*

---

## 9. Laplacian Analysis

`laplacianAnalysis` is the single validator for **both** weighted and unweighted graph Laplacians (an unweighted Laplacian is just the special case where every edge weight equals 1). It checks, in order: every entry finite, symmetry, off-diagonal entries ≤ 0, diagonal entries ≥ 0, and every row summing to zero — returning a specific rejection reason the moment any check fails, rather than silently reinterpreting an invalid matrix. If all checks pass, it extracts the edge list (see [Section 10](#10-weighted-and-unweighted-graph-data) for how weight and the weighted/unweighted flag are derived from that same loop).

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
*(`script.js:695(LAPLACIAN_TOL)`, `698-769`. Note this function calls `formatReal`, a small numeric-formatting helper defined at `script.js:300-308` — included below only because it's a direct dependency, not because formatting is itself mathematical:*

```js
function formatReal(x) {
  if (Math.abs(x) < 1e-9) x = 0;
  var rounded = Math.round(x * 1e6) / 1e6;
  if (Object.is(rounded, -0)) rounded = 0;
  var str = rounded.toFixed(6);
  str = str.replace(/(\.\d*?)0+$/, '$1').replace(/\.$/, '');
  if (str === '' || str === '-0') str = '0';
  return str;
}
```
*)*

---

## 10. Weighted and Unweighted Graph Data

This is the mathematical graph data derived from a validated Laplacian — degree, weighted degree, connected-component count, and shape classification. None of it touches pixels, SVG elements, or layout coordinates.

**Where weight actually comes from** (the mathematically meaningful part — this is the exact snippet from `laplacianAnalysis` above, repeated here because it's the crux of weighted-graph support). In plain terms: **MATHEMATICAL GRAPH DATA: `weight = |Lij|`** — that is exactly what `var w = -val;` computes below, since `val` (`= Lij`) is already known negative at this point:

```js
var val = matrix[i][j];
if (val < -LAPLACIAN_TOL) {
  var w = -val;
  edges.push([i + 1, j + 1, w]);
  if (Math.abs(w - 1) > LAPLACIAN_TOL) weighted = true;
}
```
`Lii` (the diagonal) is never read directly as a weighted degree — it is only ever *validated* against the row-sum-zero rule in Section 9. Weighted degree is instead computed independently, directly from the extracted edge list, below.

> A later, separate step (`desiredEdgeLength`, in `script.js`'s layout code, **not included in this document**) takes that same weight number and maps it to a pixel length for the on-screen edge. That mapping is visualization, not graph data — it does not change `edges`, `weight`, or any value described in this section.

**Degree and weighted degree**, computed directly from the edge list `laplacianAnalysis` returns:

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

**Connected-component count**, used by graph classification to detect a disconnected graph:

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
*(`script.js:791-802`, a union-find / disjoint-set-union over vertices.)*

**Graph classification**, a structural (degree- and edge-count-based) identification of the graph's shape — used only to produce a mathematical label (e.g. "Complete Graph K₄"); it never influences vertex/edge/weight extraction above it:

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

The "Weighted " label prefix itself (e.g. turning "Complete Graph K₄" into "Weighted Complete Graph K₄") is a one-line string concatenation done in `renderGraphSection` (`script.js:1585`, shown in Section 2) — `classifyGraph` itself is weight-agnostic and returns the same shape name either way.

---

## 11. Function Reference

| Function | Location | Purpose |
|---|---|---|
| `cx`, `cAdd`, `cSub`, `cMul`, `cDiv`, `cAbs`, `cNeg`, `cConj`, `cScaleR`, `cSqrt` | `script.js:27-47` | Complex-number arithmetic primitives used throughout the QR solver and null-space calculation. |
| `realToComplexMatrix` | `script.js:50-52` | Converts a real matrix to a complex one (imaginary part 0). |
| `toHessenberg` | `script.js:55-84` | Reduces the input matrix to upper-Hessenberg form via Householder reflections (similarity transform; preserves eigenvalues). |
| `complexGivens` | `script.js:95-104` | Computes a single complex Givens rotation for one QR step. |
| `hessenbergQRStep` | `script.js:106-132` | Applies one implicit-shift QR step to a Hessenberg matrix using Givens rotations. |
| `eigenvaluesQR` | `script.js:135-182` | Full shifted-QR eigenvalue solver with deflation; returns all eigenvalues (real or complex). |
| `internalClean` | `script.js:185-191` | Rounds a near-zero or near-integer value to a clean display/comparison value. |
| `groupEigenvalues` | `script.js:193-221` | Clusters numerically equal roots; cluster size = algebraic multiplicity. |
| `nullspace` | `script.js:224-281` | Computes a basis for the null space of `(A − λI)` via complex Gaussian elimination; basis = eigenvectors, its length = geometric multiplicity. |
| `computeEigen` | `script.js:284-295` | Main entry point: ties `eigenvaluesQR`, `groupEigenvalues`, and `nullspace` together into `{ value, multiplicity, vectors }` groups. |
| `formatReal` | `script.js:300-308` | Numeric formatting helper; included because `laplacianAnalysis` calls it when building rejection messages. |
| `laplacianAnalysis` | `script.js:708-769` | Validates a matrix as a (weighted or unweighted) graph Laplacian; extracts `edges` (with weight) and the `weighted` flag, or returns a specific rejection reason. |
| `computeDegrees` | `script.js:771-778` | Plain vertex degree (count of incident edges). |
| `computeWeightedDegrees` | `script.js:782-789` | Weighted vertex degree (sum of incident edge weights). |
| `countComponents` | `script.js:791-802` | Union-find count of connected components. |
| `toSubscript` | `script.js:805-807` | Converts a number to Unicode subscript digits, for labels like K₄. |
| `classifyGraph` | `script.js:813-842` | Structural shape classification (Kn, Pn, Cn, star, disconnected, empty, single) from edges/degrees alone. |

---

## 12. What Is NOT Part of the Mathematical Core

Intentionally excluded from this document, because it is presentation, layout, or plumbing rather than the calculation itself:

- **All HTML, CSS, colors, gradients, typography, buttons, cards, page layout, animations, and responsive styling.**
- **`formatComplex`** (`script.js:310-318`) — display formatting for eigenvalues/eigenvectors; not used by any calculation.
- **Graph *layout/positioning*** — `componentsList`, `traceChain`, `detectComponentShape`, `layoutComponent`, `computeLayoutPoints`, `desiredEdgeLength`, `forceDirectedLayout`, `edgeWeightLookup`, `layoutWeightedPath`, `layoutWeightedCycle`, `layoutWeightedStar`, `computeWeightedLayoutPoints`. These decide *where a vertex is drawn on screen* (including the weight→pixel-length mapping) — they consume the mathematical edge/weight data from Section 10 but produce visual coordinates, not graph data. The graph visualization itself is not the mathematical core.
- **SVG construction and styling** — `svgEl`, `graphVisualParams`, `edgeOpacityForCount`, `buildGraphSVG`.
- **DOM rendering and event handling** — `buildMatrixGrid`, `onCellInput`, `getCell`, `onCellKeydown`, `caretAtEnd`/`caretAtStart`, `updateLargeMatrixNotice`, `createMatrix`, `readMatrix` (DOM input parsing, not the calculation), `renderStaticMatrix`, `renderVectorBlock`, `renderResults`, `escapeHtml`, the DOM-writing bulk of `renderGraphSection`, `clearResults`, `clearMatrixValues`, `loadExample`.
- **Clipboard / file export** — `buildGraphCopyText`, `copyTextToClipboard`, `copyGraphInfo`, `downloadGraph`, `buildCopyText`, `copyResults`.
- **Top-level UI wiring** — element lookups (`document.getElementById(...)`), `showError`, and all `addEventListener` calls, except the one `runCalculation` excerpt in Section 2 kept solely to show where the two pipelines are actually invoked.
