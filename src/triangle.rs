//! Triangle and point primitives for Delaunay triangulation.

use std::fmt;

/// A 2D point with `f64` coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    /// x-coordinate.
    pub x: f64,
    /// y-coordinate.
    pub y: f64,
}

impl Point {
    /// Create a new point.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Squared distance to another point.
    pub fn dist_sq(&self, other: &Point) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    /// Distance to another point.
    pub fn distance(&self, other: &Point) -> f64 {
        self.dist_sq(other).sqrt()
    }

    /// Squared distance from this point to the circumcircle of three points.
    /// Returns the squared circumradius if the point is on the circle.
    /// Positive = outside, negative = inside.
    pub fn in_circumcircle(&self, a: &Point, b: &Point, c: &Point) -> bool {
        let ax = a.x - self.x;
        let ay = a.y - self.y;
        let bx = b.x - self.x;
        let by = b.y - self.y;
        let cx = c.x - self.x;
        let cy = c.y - self.y;

        let det = ax * (by * (cx * cx + cy * cy) - cy * (bx * bx + by * by))
            - ay * (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by))
            + (ax * ax + ay * ay) * (bx * cy - by * cx);

        // For CCW triangles, positive determinant means inside
        let orient = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        if orient > 0.0 { det > 0.0 } else { det < 0.0 }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.4}, {:.4})", self.x, self.y)
    }
}

/// A triangle defined by three vertex indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Triangle {
    /// First vertex index.
    pub a: usize,
    /// Second vertex index.
    pub b: usize,
    /// Third vertex index.
    pub c: usize,
}

impl Triangle {
    /// Create a new triangle from three vertex indices.
    pub fn new(a: usize, b: usize, c: usize) -> Self {
        Self { a, b, c }
    }

    /// Check if the triangle contains a given vertex index.
    pub fn contains_vertex(&self, v: usize) -> bool {
        self.a == v || self.b == v || self.c == v
    }

    /// Get the three vertex indices.
    pub fn vertices(&self) -> [usize; 3] {
        [self.a, self.b, self.c]
    }

    /// Does this triangle share an edge with another?
    pub fn shares_edge(&self, other: &Triangle) -> bool {
        let sv = self.vertices();
        let ov = other.vertices();
        let shared = sv.iter().filter(|v| ov.contains(v)).count();
        shared == 2
    }

    /// Compute circumcircle center and squared radius for the triangle.
    pub fn circumcircle(&self, points: &[Point]) -> (Point, f64) {
        let pa = &points[self.a];
        let pb = &points[self.b];
        let pc = &points[self.c];

        let d = 2.0 * (pa.x * (pb.y - pc.y) + pb.x * (pc.y - pa.y) + pc.x * (pa.y - pb.y));

        if d.abs() < 1e-12 {
            // Degenerate triangle
            return (Point::new(0.0, 0.0), f64::MAX);
        }

        let ux = ((pa.x * pa.x + pa.y * pa.y) * (pb.y - pc.y)
            + (pb.x * pb.x + pb.y * pb.y) * (pc.y - pa.y)
            + (pc.x * pc.x + pc.y * pc.y) * (pa.y - pb.y))
            / d;

        let uy = ((pa.x * pa.x + pa.y * pa.y) * (pc.x - pb.x)
            + (pb.x * pb.x + pb.y * pb.y) * (pa.x - pc.x)
            + (pc.x * pc.x + pc.y * pc.y) * (pb.x - pa.x))
            / d;

        let center = Point::new(ux, uy);
        let r_sq = center.dist_sq(pa);
        (center, r_sq)
    }

    /// Compute signed area (positive if CCW).
    pub fn signed_area(&self, points: &[Point]) -> f64 {
        let pa = &points[self.a];
        let pb = &points[self.b];
        let pc = &points[self.c];
        0.5 * ((pb.x - pa.x) * (pc.y - pa.y) - (pc.x - pa.x) * (pb.y - pa.y))
    }

    /// Area (always positive).
    pub fn area(&self, points: &[Point]) -> f64 {
        self.signed_area(points).abs()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_new() {
        let p = Point::new(1.0, 2.0);
        assert_eq!(p.x, 1.0);
        assert_eq!(p.y, 2.0);
    }

    #[test]
    fn test_distance() {
        let a = Point::new(0.0, 0.0);
        let b = Point::new(3.0, 4.0);
        assert!((a.distance(&b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_triangle_new() {
        let t = Triangle::new(0, 1, 2);
        assert_eq!(t.vertices(), [0, 1, 2]);
    }

    #[test]
    fn test_contains_vertex() {
        let t = Triangle::new(0, 1, 2);
        assert!(t.contains_vertex(1));
        assert!(!t.contains_vertex(3));
    }

    #[test]
    fn test_shares_edge() {
        let t1 = Triangle::new(0, 1, 2);
        let t2 = Triangle::new(1, 2, 3);
        let t3 = Triangle::new(3, 4, 5);
        assert!(t1.shares_edge(&t2));
        assert!(!t1.shares_edge(&t3));
    }

    #[test]
    fn test_circumcircle() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(0.5, 0.866025),
        ];
        let t = Triangle::new(0, 1, 2);
        let (center, r_sq) = t.circumcircle(&pts);
        assert!(r_sq > 0.0);
        // All vertices should be equidistant from center
        let d0 = center.dist_sq(&pts[0]);
        let d1 = center.dist_sq(&pts[1]);
        let d2 = center.dist_sq(&pts[2]);
        assert!((d0 - r_sq).abs() < 1e-6);
        assert!((d1 - r_sq).abs() < 1e-6);
        assert!((d2 - r_sq).abs() < 1e-6);
    }

    #[test]
    fn test_area() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(4.0, 0.0),
            Point::new(0.0, 3.0),
        ];
        let t = Triangle::new(0, 1, 2);
        assert!((t.area(&pts) - 6.0).abs() < 1e-10);
    }
}
