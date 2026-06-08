# delaunay-triang-rs

[![crates.io](https://img.shields.io/crates/v/delaunay-triang-rs.svg)](https://crates.io/crates/delaunay-triang-rs)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Delaunay triangulation and Voronoi diagrams in pure Rust — Bowyer-Watson incremental insertion, edge flip, and Guibas-Stolfi quad-edge.

## The Problem

Given a set of 2D points, you need the "best" triangulation — one that avoids skinny triangles and is uniquely defined by the input points. This is the Delaunay triangulation: for every triangle, no other input point lies inside its circumcircle.

The dual of the Delaunay triangulation is the **Voronoi diagram**, which partitions the plane into regions closest to each site. Both structures are fundamental in computational geometry: mesh generation, nearest-neighbor queries, terrain modeling, and spatial interpolation.

## The Insight

The Delaunay property has a clean geometric characterization: a triangulation is Delaunay if and only if for every edge, the sum of the opposite angles in the two adjacent triangles is ≤ 180°. Equivalently, no point lies inside any triangle's circumcircle.

This property enables two algorithmic approaches:

1. **Incremental insertion (Bowyer-Watson):** Add points one at a time. For each new point, find all triangles whose circumcircle contains it (the "cavity"), delete them, and re-triangulate the resulting polygonal hole with the new point.
2. **Edge flipping:** Start with any triangulation. Find edges that violate the Delaunay property (opposite angles sum > 180°) and flip them. Repeat until no more flips are needed. Guaranteed to terminate because each flip increases the minimum angle.

## How It Works

### Bowyer-Watson (`bowyer_watson`)

1. Create a **super-triangle** large enough to contain all input points.
2. For each point:
   - Find all triangles whose circumcircle contains the point.
   - Identify the boundary polygon of the cavity (edges shared with exactly one bad triangle).
   - Delete bad triangles; create new triangles from each boundary edge to the new point.
3. Remove all triangles that reference super-triangle vertices.

The super-triangle trick ensures the point always lies inside some triangle's circumcircle, so the cavity is never empty.

### Edge Flip (`edge_flip`)

Given an existing triangulation:
1. For each pair of adjacent triangles sharing an edge (a, b) with opposite vertices c and d.
2. If point d lies inside the circumcircle of triangle (a, b, c), the edge is **illegal**.
3. Flip: replace edge (a, b) with edge (c, d), creating triangles (c, d, a) and (c, d, b).
4. Repeat until no illegal edges remain.

Bounded by O(n²) flips in the worst case, but typically O(n) for well-behaved inputs.

### Quad-Edge (`quad_edge`)

The Guibas-Stolfi data structure represents a planar subdivision using **quad-edges** — each logical edge is stored as four directed half-edges (original, rotated 90°, symmetric, rotated 270°). Operations:

- `onext(e)`: Next edge CCW around origin
- `sym(e)`: Opposite direction
- `rot(e)` / `rot_inv(e)`: Rotate 90° CCW / CW
- `splice(a, b)`: Topological operation that connects/disconnects edge rings

Splice is the fundamental primitive — `make_edge` and `delete_edge` are built from it.

### Voronoi (`voronoi`)

Computed from the Delaunay triangulation dual:
1. Compute the circumcenter of each Delaunay triangle.
2. For each pair of adjacent triangles sharing an edge, connect their circumcenters with a Voronoi edge.
3. The shared edge's endpoints identify which two sites the Voronoi edge separates.

## Usage

```rust
use delaunay_triang::{Point, Triangle};
use delaunay_triang::bowyer_watson::triangulate;
use delaunay_triang::voronoi::VoronoiDiagram;
use delaunay_triang::edge_flip::edge_flip_delaunay;
use delaunay_triang::quad_edge::QuadEdgeSubdivision;

// --- Bowyer-Watson triangulation ---
let points = vec![
    Point::new(0.0, 0.0),
    Point::new(1.0, 0.0),
    Point::new(0.5, 1.0),
    Point::new(1.5, 1.0),
];
let triangles = triangulate(&points);
assert_eq!(triangles.len(), 2);

// --- Voronoi diagram from Delaunay ---
let voronoi = VoronoiDiagram::from_delaunay(&points, &triangles);
for edge in &voronoi.edges {
    println!("Voronoi edge: {} -> {} (sites: {}, {})",
             edge.from, edge.to, edge.site_left, edge.site_right);
}
// Find edges for a specific site
let site_edges = voronoi.edges_for_site(0);

// --- Edge flip ---
let points = vec![
    Point::new(0.0, 0.0),
    Point::new(1.0, 0.0),
    Point::new(1.0, 1.0),
    Point::new(0.0, 1.0),
];
let mut tris = vec![
    Triangle::new(0, 1, 2),
    Triangle::new(0, 2, 3),
];
edge_flip_delaunay(&points, &mut tris);

// --- Quad-edge subdivision ---
let mut qes = QuadEdgeSubdivision::new();
let a = qes.add_point(Point::new(0.0, 0.0));
let b = qes.add_point(Point::new(1.0, 0.0));
let e = qes.make_edge(a, b);
assert_eq!(qes.origin(e), Some(a));
assert_eq!(qes.dest(e), Some(b));
```

## Module Map

| Module | Description |
|---|---|
| `triangle` | `Point` (2D f64) and `Triangle` (three vertex indices). Circumcircle, area, edge sharing. |
| `bowyer_watson` | Incremental Delaunay triangulation via super-triangle approach |
| `edge_flip` | Convert any triangulation to Delaunay via iterative edge flipping |
| `quad_edge` | Guibas-Stolfi quad-edge planar subdivision data structure |
| `voronoi` | Voronoi diagram construction from Delaunay dual |

Re-exports: `Point` and `Triangle` from the crate root.

## Design Decisions

- **Super-triangle (not bounding box).** The Bowyer-Watson implementation uses a large triangle (20× the bounding box dimension) rather than an infinite sentinel. This is simpler but requires post-processing to remove super-triangle vertices.
- **Brute-force cavity search.** For each new point, all triangles are checked against the circumcircle. This is O(n) per insertion, yielding O(n²) total. A spatial index (e.g., a history DAG or point location structure) would reduce this to O(n log n), but adds significant complexity.
- **`f64` everywhere.** No exact arithmetic. The circumcircle test uses a tolerance (`1e-10`) to handle near-degenerate cases. This means the Delaunay property is approximately satisfied for nearly-cocircular points.
- **Vertex indices, not references.** Triangles store `usize` indices into the point array. This avoids lifetime issues and makes it easy to serialize, but requires the caller to maintain the point array alongside the triangles.
- **Quad-edge as a standalone structure.** The quad-edge module doesn't implement Delaunay construction — it's a general-purpose planar subdivision data structure that could be used to implement incremental Delaunay, Voronoi, or other planar algorithms.

## Complexity

| Algorithm | Time | Space |
|---|---|---|
| Bowyer-Watson | O(n²) worst case, O(n log n) typical with point location | O(n) |
| Edge flip | O(n²) worst case | O(n) in-place |
| Voronoi from Delaunay | O(t²) where t = triangles | O(t) |
| Quad-edge operations | O(1) per `splice`, `make_edge` | O(e) where e = edges |

## Limitations

- **No 3D.** All algorithms are strictly 2D.
- **Numerical robustness.** Uses `f64` with epsilon comparisons. Nearly-cocircular or collinear points may produce incorrect triangulations. For production geometry software, use exact arithmetic (e.g., `rust-geo`'s robust predicates).
- **No constrained Delaunay.** Can't enforce required edges (e.g., polygon boundaries). The triangulation is always unconstrained.
- **Quadratic Bowyer-Watson.** The brute-force cavity search makes this unsuitable for > 10,000 points without modification. For large point sets, add a point-location structure.
- **Voronoi doesn't clip to bounding box.** Voronoi edges extend to circumcenters of boundary triangles, which may lie outside the convex hull. Clipping requires additional ray-boundary intersection logic.
- **No parallelism.** All algorithms are single-threaded.

## Status

Published to [crates.io](https://crates.io/crates/delaunay-triang-rs). Covers the fundamental Delaunay/Voronoi algorithms with clean module separation. Suitable for educational use, small-scale mesh generation, and geometry prototyping. For production computational geometry, consider `spade` (Rust) or `CGAL` (C++) for robust predicates and optimized algorithms.
