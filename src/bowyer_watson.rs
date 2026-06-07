//! Bowyer-Watson algorithm for Delaunay triangulation.
//!
//! Incrementally inserts points into a triangulation, maintaining the
//! Delaunay property by removing and re-triangulating affected triangles.

use crate::triangle::{Point, Triangle};

/// Compute the Delaunay triangulation of a set of points using the Bowyer-Watson algorithm.
///
/// Returns a list of triangles as vertex index triples.
/// Uses a super-triangle approach for initialization.
pub fn triangulate(points: &[Point]) -> Vec<Triangle> {
    if points.len() < 3 {
        return Vec::new();
    }

    // Create super-triangle that encompasses all points
    let (min_x, max_x, min_y, max_y) = bounding_box(points);
    let dx = max_x - min_x;
    let dy = max_y - min_y;
    let d = dx.max(dy).max(1.0);
    let cx = (min_x + max_x) / 2.0;
    let cy = (min_y + max_y) / 2.0;

    // Super-triangle vertices (at indices n, n+1, n+2)
    let st_a = Point::new(cx - 20.0 * d, cy - d);
    let st_b = Point::new(cx, cy + 20.0 * d);
    let st_c = Point::new(cx + 20.0 * d, cy - d);

    let mut all_points = points.to_vec();
    let n = all_points.len();
    all_points.push(st_a);
    all_points.push(st_b);
    all_points.push(st_c);

    let mut triangles = vec![Triangle::new(n, n + 1, n + 2)];

    // Insert each point
    for i in 0..n {
        let p = all_points[i];
        let mut bad_triangles = Vec::new();

        // Find all triangles whose circumcircle contains the new point
        for (ti, t) in triangles.iter().enumerate() {
            let (center, r_sq) = t.circumcircle(&all_points);
            if center.dist_sq(&p) < r_sq + 1e-10 {
                bad_triangles.push(ti);
            }
        }

        // Find boundary of the polygonal hole
        let mut polygon: Vec<(usize, usize)> = Vec::new();
        for &ti in &bad_triangles {
            let t = triangles[ti];
            let edges = [(t.a, t.b), (t.b, t.c), (t.c, t.a)];
            for edge in edges {
                let mut shared = false;
                'inner: for &tj in &bad_triangles {
                    if tj == ti { continue; }
                    let ot = triangles[tj];
                    let oedges = [(ot.a, ot.b), (ot.b, ot.c), (ot.c, ot.a)];
                    for oe in &oedges {
                        if (edge.0 == oe.0 && edge.1 == oe.1) || (edge.0 == oe.1 && edge.1 == oe.0) {
                            shared = true;
                            break 'inner;
                        }
                    }
                }
                if !shared {
                    polygon.push(edge);
                }
            }
        }

        // Remove bad triangles (in reverse order to maintain indices)
        let mut bad_sorted = bad_triangles;
        bad_sorted.sort_unstable();
        for &ti in bad_sorted.iter().rev() {
            triangles.remove(ti);
        }

        // Create new triangles from polygon edges to the new point
        for edge in &polygon {
            triangles.push(Triangle::new(edge.0, edge.1, i));
        }
    }

    // Remove triangles that share vertices with the super-triangle
    triangles.retain(|t| {
        t.a < n && t.b < n && t.c < n
    });

    triangles
}

/// Compute the bounding box of a set of points.
fn bounding_box(points: &[Point]) -> (f64, f64, f64, f64) {
    let mut min_x = f64::MAX;
    let mut max_x = f64::MIN;
    let mut min_y = f64::MAX;
    let mut max_y = f64::MIN;
    for p in points {
        min_x = min_x.min(p.x);
        max_x = max_x.max(p.x);
        min_y = min_y.min(p.y);
        max_y = max_y.max(p.y);
    }
    (min_x, max_x, min_y, max_y)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_triangulation() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(0.5, 1.0),
            Point::new(1.5, 1.0),
        ];
        let tris = triangulate(&pts);
        assert!(tris.len() >= 2);
        // Euler's formula: T = 2n - 2 - h for n points with h on convex hull
        // For 4 points with convex hull of 4: T = 2
        assert_eq!(tris.len(), 2);
    }

    #[test]
    fn test_three_points() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(0.0, 1.0),
        ];
        let tris = triangulate(&pts);
        assert_eq!(tris.len(), 1);
    }

    #[test]
    fn test_two_points() {
        let pts = vec![Point::new(0.0, 0.0), Point::new(1.0, 1.0)];
        let tris = triangulate(&pts);
        assert!(tris.is_empty());
    }

    #[test]
    fn test_collinear_points() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(2.0, 0.0),
            Point::new(3.0, 0.0),
        ];
        let tris = triangulate(&pts);
        // Collinear points may produce degenerate triangles or none
        for t in &tris {
            assert!(t.area(&pts) >= -1e-6);
        }
    }

    #[test]
    fn test_square_points() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(1.0, 1.0),
            Point::new(0.0, 1.0),
        ];
        let tris = triangulate(&pts);
        assert_eq!(tris.len(), 2);
        // Total area should be ~1.0
        let total_area: f64 = tris.iter().map(|t| t.area(&pts)).sum();
        assert!((total_area - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_delaunay_property() {
        // Verify no point is inside any circumcircle
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(2.0, 0.0),
            Point::new(1.0, 2.0),
            Point::new(0.5, 0.5),
            Point::new(1.5, 0.5),
        ];
        let tris = triangulate(&pts);
        for t in &tris {
            let (center, r_sq) = t.circumcircle(&pts);
            for (i, p) in pts.iter().enumerate() {
                if t.contains_vertex(i) { continue; }
                let d = center.dist_sq(p);
                assert!(d >= r_sq - 1e-6, "Point {} inside circumcircle of triangle {:?}", i, t);
            }
        }
    }

    #[test]
    fn test_all_vertices_valid() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(0.5, 1.0),
            Point::new(-0.5, 0.5),
            Point::new(1.5, 0.5),
        ];
        let tris = triangulate(&pts);
        let n = pts.len();
        for t in &tris {
            assert!(t.a < n);
            assert!(t.b < n);
            assert!(t.c < n);
        }
    }

    #[test]
    fn test_grid_triangulation() {
        let mut pts = Vec::new();
        for i in 0..4 {
            for j in 0..4 {
                pts.push(Point::new(i as f64, j as f64));
            }
        }
        let tris = triangulate(&pts);
        // 4x4 grid = 16 points, expected ~18 triangles
        assert!(tris.len() >= 10);
    }
}
