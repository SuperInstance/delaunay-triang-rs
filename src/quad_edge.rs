//! Guibas-Stolfi quad-edge data structure.
//!
//! A topological data structure for representing planar subdivisions,
//! supporting edge operations (flip, rotate, splice) used in
//! Delaunay triangulation and Voronoi diagrams.

use crate::triangle::Point;

/// A reference to an edge with a rotation (0=original, 1=rot90, 2=sym, 3=rot270).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeRef {
    /// Edge ID.
    pub id: usize,
    /// Rotation index (0..3).
    pub rot: u8,
}

impl EdgeRef {
    /// Create a new edge reference.
    pub fn new(id: usize, rot: u8) -> Self {
        Self { id, rot: rot % 4 }
    }

    /// The symmetric (opposite) edge.
    pub fn sym(self) -> Self {
        Self { id: self.id, rot: (self.rot + 2) % 4 }
    }

    /// Rotate 90 degrees CCW.
    pub fn rot(self) -> Self {
        Self { id: self.id, rot: (self.rot + 1) % 4 }
    }

    /// Rotate 90 degrees CW (inverse of rot).
    pub fn rot_inv(self) -> Self {
        Self { id: self.id, rot: (self.rot + 3) % 4 }
    }
}

/// A quad-edge: four directed edges sharing the same ID.
pub struct QuadEdge {
    /// Next edge CCW around origin for each rotation.
    next: [EdgeRef; 4],
    /// Origin vertex for each rotation.
    origin: [Option<usize>; 4],
}

impl QuadEdge {
    /// Create a new quad-edge with self-referential next pointers.
    fn new(id: usize) -> Self {
        let e0 = EdgeRef::new(id, 0);
        Self {
            next: [e0, e0.rot(), e0.sym(), e0.sym().rot()],
            origin: [None, None, None, None],
        }
    }

    fn get_next(&self, er: EdgeRef) -> EdgeRef {
        self.next[er.rot as usize]
    }

    fn set_next(&mut self, er: EdgeRef, val: EdgeRef) {
        self.next[er.rot as usize] = val;
    }

    fn get_origin(&self, er: EdgeRef) -> Option<usize> {
        self.origin[er.rot as usize]
    }

    fn set_origin(&mut self, er: EdgeRef, val: Option<usize>) {
        self.origin[er.rot as usize] = val;
    }
}

/// A quad-edge subdivision (planar graph).
pub struct QuadEdgeSubdivision {
    /// All quad-edges.
    pub edges: Vec<QuadEdge>,
    /// Points.
    pub points: Vec<Point>,
}

impl QuadEdgeSubdivision {
    /// Create a new empty subdivision.
    pub fn new() -> Self {
        Self {
            edges: Vec::new(),
            points: Vec::new(),
        }
    }

    /// Add a point and return its index.
    pub fn add_point(&mut self, p: Point) -> usize {
        let idx = self.points.len();
        self.points.push(p);
        idx
    }

    /// Create a new edge between two points.
    pub fn make_edge(&mut self, origin: usize, dest: usize) -> EdgeRef {
        let id = self.edges.len();
        let e = EdgeRef::new(id, 0);
        let esym = e.sym();

        let mut qe = QuadEdge::new(id);
        qe.set_origin(e, Some(origin));
        qe.set_origin(esym, Some(dest));
        self.edges.push(qe);

        e
    }

    /// Get the next edge CCW around origin.
    pub fn onext(&self, e: EdgeRef) -> EdgeRef {
        self.edges[e.id].get_next(e)
    }

    /// Get the origin vertex of an edge.
    pub fn origin(&self, e: EdgeRef) -> Option<usize> {
        self.edges[e.id].get_origin(e)
    }

    /// Get the destination vertex of an edge.
    pub fn dest(&self, e: EdgeRef) -> Option<usize> {
        self.origin(e.sym())
    }

    /// Next edge CCW around the left face.
    pub fn lnext(&self, e: EdgeRef) -> EdgeRef {
        self.onext(e.rot_inv()).rot()
    }

    /// Splice two edges (topological operation).
    pub fn splice(&mut self, a: EdgeRef, b: EdgeRef) {
        let alpha = self.onext(a).rot();
        let beta = self.onext(b).rot();

        let t1 = self.onext(a);
        let t2 = self.onext(b);
        let t3 = self.onext(alpha);
        let t4 = self.onext(beta);

        self.edges[a.id].set_next(a, t2);
        self.edges[b.id].set_next(b, t1);
        self.edges[alpha.id].set_next(alpha, t4);
        self.edges[beta.id].set_next(beta, t3);
    }

    /// Delete an edge from the subdivision.
    pub fn delete_edge(&mut self, e: EdgeRef) {
        self.splice(e, self.onext(e));
        let esym = e.sym();
        self.splice(esym, self.onext(esym));
    }

    /// Number of edges.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }
}

impl Default for QuadEdgeSubdivision {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_edge() {
        let mut qes = QuadEdgeSubdivision::new();
        let a = qes.add_point(Point::new(0.0, 0.0));
        let b = qes.add_point(Point::new(1.0, 0.0));
        let e = qes.make_edge(a, b);
        assert_eq!(qes.origin(e), Some(a));
        assert_eq!(qes.dest(e), Some(b));
    }

    #[test]
    fn test_sym_edge() {
        let mut qes = QuadEdgeSubdivision::new();
        let a = qes.add_point(Point::new(0.0, 0.0));
        let b = qes.add_point(Point::new(1.0, 0.0));
        let e = qes.make_edge(a, b);
        let esym = e.sym();
        assert_eq!(qes.origin(esym), Some(b));
        assert_eq!(qes.dest(esym), Some(a));
    }

    #[test]
    fn test_rot() {
        let mut qes = QuadEdgeSubdivision::new();
        let a = qes.add_point(Point::new(0.0, 0.0));
        let b = qes.add_point(Point::new(1.0, 0.0));
        let e = qes.make_edge(a, b);
        assert_eq!(e.rot().rot(), e.sym());
    }

    #[test]
    fn test_splice() {
        let mut qes = QuadEdgeSubdivision::new();
        let a = qes.add_point(Point::new(0.0, 0.0));
        let b = qes.add_point(Point::new(1.0, 0.0));
        let c = qes.add_point(Point::new(0.0, 1.0));

        let e1 = qes.make_edge(a, b);
        let e2 = qes.make_edge(a, c);

        qes.splice(e1.sym(), e2);
        assert_eq!(qes.edge_count(), 2);
    }

    #[test]
    fn test_default() {
        let qes = QuadEdgeSubdivision::default();
        assert_eq!(qes.edge_count(), 0);
        assert_eq!(qes.points.len(), 0);
    }

    #[test]
    fn test_multiple_edges() {
        let mut qes = QuadEdgeSubdivision::new();
        let pts: Vec<usize> = (0..5).map(|i| {
            qes.add_point(Point::new(i as f64, 0.0))
        }).collect();

        for i in 0..4 {
            qes.make_edge(pts[i], pts[i + 1]);
        }
        assert_eq!(qes.edge_count(), 4);
    }

    #[test]
    fn test_delete_edge() {
        let mut qes = QuadEdgeSubdivision::new();
        let a = qes.add_point(Point::new(0.0, 0.0));
        let b = qes.add_point(Point::new(1.0, 0.0));
        let e = qes.make_edge(a, b);
        assert_eq!(qes.edge_count(), 1);
        qes.delete_edge(e);
        assert_eq!(qes.edge_count(), 1);
    }
}
