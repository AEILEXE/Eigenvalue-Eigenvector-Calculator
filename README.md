# Eigenvalue & Eigenvector Calculator

A browser-based calculator that computes the eigenvalues and eigenvectors of any square matrix. Enter a matrix size, fill in the values, and get the eigenvalues, their multiplicities, and their corresponding eigenvectors — all calculated instantly in your browser, with no data ever sent to a server.

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

In linear algebra, an eigenvalue and eigenvector pair describes a direction that a matrix only stretches or shrinks, without changing — a property used throughout engineering, physics, computer graphics, and data science (for example, in vibration analysis, stability analysis, and principal component analysis). This calculator handles that computation for you, including matrices whose eigenvalues turn out to be complex numbers.

---

## Features

- Square matrix input for any size from 1×1 to 100×100
- Eigenvalue calculation, including complex (non-real) eigenvalues
- Eigenvector calculation for every eigenvalue
- Algebraic multiplicity (how many times an eigenvalue repeats as a root)
- Geometric multiplicity (how many independent eigenvectors an eigenvalue actually has)
- A "Defective" indicator when an eigenvalue's algebraic and geometric multiplicities differ
- A built-in example matrix you can load with one click
- One-click clearing of the current matrix
- One-click copying of the full results as plain text
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

## Interface Guide

| Control | Purpose |
|---|---|
| Matrix Size (n) | Sets the dimension of the square matrix (1–100) |
| Create Matrix | Builds the matrix input grid for the chosen size |
| Calculate Eigenvalues & Eigenvectors | Computes and displays the eigenvalues and eigenvectors |
| Load Example | Loads a built-in 2×2 example matrix |
| Clear Matrix | Empties all cell values in the current matrix |
| Copy Results | Copies the matrix and calculated results to the clipboard |

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
- **JavaScript (vanilla, no frameworks)** — all interactivity and calculations

No external math or linear algebra library is used. The eigenvalue and eigenvector computations (Hessenberg reduction, the shifted QR algorithm with complex-number support, and null-space calculation) are implemented directly in the project's own JavaScript. There is no backend, server, or database — the entire calculator runs client-side in the browser.

---

## How It Works

1. The user selects a matrix size and clicks Create Matrix (or presses Enter).
2. The website generates an empty grid of input cells for that size.
3. The user enters a numeric value into every cell.
4. When Calculate is clicked, JavaScript reads and validates the matrix values.
5. The matrix is reduced to a simpler equivalent form, then processed with a numerical algorithm to find its eigenvalues.
6. For each eigenvalue found, the calculator determines its multiplicity and computes its associated eigenvector(s).
7. The eigenvalues, multiplicities, and eigenvectors are formatted and displayed on the page.

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
