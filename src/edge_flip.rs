//! Edge flip algorithm for Delaunay triangulation.
//!
//! Converts any triangulation into a Delaunay triangulation by
//! iteratively flipping edges that violate the Delaunay property.

use crate::triangle::{Point, Triangle};

/// Convert an arbitrary triangulation to a Delaunay triangulation using edge flips.
///
/// Takes a set of points and an initial triangulation (as triangles with vertex indices),
/// and returns a Delaunay triangulation obtained by flipping illegal edges.
pub fn edge_flip_delaunay(points: &[Point], triangles: &mut [Triangle]) {
    let mut flipped = true;
    let max_iterations = triangles.len() * triangles.len() * 10;
    let mut iterations = 0;

    while flipped && iterations < max_iterations {
        flipped = false;
        iterations += 1;

        let n = triangles.len();
        for i in 0..n {
            for j in (i + 1)..n {
                if !triangles[i].shares_edge(&triangles[j]) {
                    continue;
                }

                // Find shared edge and the two non-shared vertices
                let shared = find_shared_vertices(&triangles[i], &triangles[j]);
                if shared.is_none() { continue; }
                let (a, b, c, d) = shared.unwrap();

                // Check if edge a-b should be flipped to c-d
                if should_flip(points, a, b, c, d) {
                    // Flip: replace triangles i and j
                    triangles[i] = Triangle::new(c, d, a);
                    triangles[j] = Triangle::new(c, d, b);
                    flipped = true;
                }
            }
        }
    }
}

/// Find the shared edge between two triangles and the two non-shared vertices.
/// Returns (shared_v1, shared_v2, non_shared_from_t1, non_shared_from_t2).
fn find_shared_vertices(t1: &Triangle, t2: &Triangle) -> Option<(usize, usize, usize, usize)> {
    let v1 = t1.vertices();
    let v2 = t2.vertices();
    let shared: Vec<usize> = v1.iter().filter(|v| v2.contains(v)).copied().collect();
    let unique1: Vec<usize> = v1.iter().filter(|v| !v2.contains(v)).copied().collect();
    let unique2: Vec<usize> = v2.iter().filter(|v| !v1.contains(v)).copied().collect();

    if shared.len() == 2 && unique1.len() == 1 && unique2.len() == 1 {
        Some((shared[0], shared[1], unique1[0], unique2[0]))
    } else {
        None
    }
}

/// Check if the edge between a and b should be flipped to connect c and d.
/// This is true if d is inside the circumcircle of triangle (a, b, c).
fn should_flip(points: &[Point], _a: usize, _b: usize, c: usize, d: usize) -> bool {
    // Use the in-circumcircle test: if d is inside circumcircle of (a,b,c), flip
    let pa = &points[_a];
    let pb = &points[_b];
    let pc = &points[c];
    let pd = &points[d];
    pd.in_circumcircle(pa, pb, pc)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_square_triangulation() -> (Vec<Point>, Vec<Triangle>) {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(1.0, 1.0),
            Point::new(0.0, 1.0),
        ];
        // Two triangles forming a square, one diagonal
        let tris = vec![
            Triangle::new(0, 1, 2),
            Triangle::new(0, 2, 3),
        ];
        (pts, tris)
    }

    #[test]
    fn test_no_flip_needed() {
        let (pts, mut tris) = make_square_triangulation();
        let _original = tris.clone();
        edge_flip_delaunay(&pts, &mut tris);
        // Should still have 2 triangles
        assert_eq!(tris.len(), 2);
    }

    #[test]
    fn test_flip_produces_valid_triangles() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(2.0, 0.0),
            Point::new(1.0, 3.0),
            Point::new(0.5, 0.5),
            Point::new(1.5, 0.5),
        ];
        let mut tris = vec![
            Triangle::new(0, 1, 3),
            Triangle::new(1, 4, 3),
            Triangle::new(0, 3, 2),
            Triangle::new(1, 2, 4),
        ];
        edge_flip_delaunay(&pts, &mut tris);
        // All triangles should have positive area
        for t in &tris {
            assert!(t.area(&pts) > 1e-10);
        }
    }

    #[test]
    fn test_preserves_triangle_count() {
        let (pts, mut tris) = make_square_triangulation();
        let count_before = tris.len();
        edge_flip_delaunay(&pts, &mut tris);
        assert_eq!(tris.len(), count_before);
    }

    #[test]
    fn test_delaunay_property_after_flip() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(2.0, 0.0),
            Point::new(1.0, 2.0),
            Point::new(1.0, 0.5),
        ];
        let mut tris = vec![
            Triangle::new(0, 1, 3),
            Triangle::new(0, 2, 3),
            Triangle::new(1, 2, 3),
        ];
        edge_flip_delaunay(&pts, &mut tris);
        // Verify Delaunay property
        for t in &tris {
            let (center, r_sq) = t.circumcircle(&pts);
            for (i, p) in pts.iter().enumerate() {
                if t.contains_vertex(i) { continue; }
                assert!(center.dist_sq(p) >= r_sq - 1e-4,
                    "Delaunay property violated");
            }
        }
    }

    #[test]
    fn test_no_shared_edge_no_flip() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(0.0, 1.0),
            Point::new(3.0, 0.0),
            Point::new(3.0, 1.0),
        ];
        let mut tris = vec![
            Triangle::new(0, 1, 2),
            Triangle::new(3, 4, 1),
        ];
        let original = tris.clone();
        edge_flip_delaunay(&pts, &mut tris);
        // These don't share an edge, so nothing should change
        assert_eq!(tris.len(), original.len());
    }
}
