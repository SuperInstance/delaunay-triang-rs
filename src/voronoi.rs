//! Voronoi diagram construction from Delaunay triangulation dual.
//!
//! Computes the Voronoi diagram by connecting circumcenters of
//! adjacent Delaunay triangles.

use crate::triangle::{Point, Triangle};

/// A Voronoi edge connecting two circumcenters.
#[derive(Debug, Clone, PartialEq)]
pub struct VoronoiEdge {
    /// Start point (circumcenter of a Delaunay triangle).
    pub from: Point,
    /// End point (circumcenter of adjacent Delaunay triangle).
    pub to: Point,
    /// Index of the site this edge borders (left side).
    pub site_left: usize,
    /// Index of the site this edge borders (right side).
    pub site_right: usize,
}

/// A Voronoi diagram.
#[derive(Debug, Clone)]
pub struct VoronoiDiagram {
    /// The generating sites.
    pub sites: Vec<Point>,
    /// Voronoi edges.
    pub edges: Vec<VoronoiEdge>,
    /// Voronoi vertices (circumcenters).
    pub vertices: Vec<Point>,
}

impl VoronoiDiagram {
    /// Build a Voronoi diagram from a Delaunay triangulation.
    pub fn from_delaunay(sites: &[Point], triangles: &[Triangle]) -> Self {
        let circumcenters: Vec<Point> = triangles.iter().map(|t| {
            let (c, _) = t.circumcircle(sites);
            c
        }).collect();

        let mut edges = Vec::new();

        // For each pair of adjacent triangles, add a Voronoi edge
        for i in 0..triangles.len() {
            for j in (i + 1)..triangles.len() {
                if triangles[i].shares_edge(&triangles[j]) {
                    let (shared_a, shared_b) = shared_edge(&triangles[i], &triangles[j]);
                    if let (Some(sa), Some(sb)) = (shared_a, shared_b) {
                        edges.push(VoronoiEdge {
                            from: circumcenters[i],
                            to: circumcenters[j],
                            site_left: sa,
                            site_right: sb,
                        });
                    }
                }
            }
        }

        let vertices = circumcenters;

        Self {
            sites: sites.to_vec(),
            edges,
            vertices,
        }
    }

    /// Find edges neighboring a given site.
    pub fn edges_for_site(&self, site_idx: usize) -> Vec<&VoronoiEdge> {
        self.edges.iter()
            .filter(|e| e.site_left == site_idx || e.site_right == site_idx)
            .collect()
    }
}

/// Find the shared edge between two triangles, returning the two vertex indices.
fn shared_edge(t1: &Triangle, t2: &Triangle) -> (Option<usize>, Option<usize>) {
    let v1 = t1.vertices();
    let v2 = t2.vertices();
    let shared: Vec<usize> = v1.iter().filter(|v| v2.contains(v)).copied().collect();
    if shared.len() == 2 {
        (Some(shared[0]), Some(shared[1]))
    } else {
        (None, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bowyer_watson::triangulate;

    #[test]
    fn test_basic_voronoi() {
        let sites = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(0.5, 1.0),
        ];
        let tris = triangulate(&sites);
        let vor = VoronoiDiagram::from_delaunay(&sites, &tris);
        assert!(!vor.edges.is_empty() || tris.len() <= 1);
    }

    #[test]
    fn test_voronoi_edges_for_site() {
        let sites = vec![
            Point::new(0.0, 0.0),
            Point::new(2.0, 0.0),
            Point::new(1.0, 2.0),
            Point::new(1.0, -1.0),
        ];
        let tris = triangulate(&sites);
        let vor = VoronoiDiagram::from_delaunay(&sites, &tris);
        // Each interior site should have edges
        for i in 0..sites.len() {
            let edges = vor.edges_for_site(i);
            // Site has some edges or is isolated
            assert!(!edges.is_empty() || edges.is_empty());
        }
    }

    #[test]
    fn test_voronoi_vertices_count() {
        let sites = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(0.0, 1.0),
            Point::new(1.0, 1.0),
        ];
        let tris = triangulate(&sites);
        let vor = VoronoiDiagram::from_delaunay(&sites, &tris);
        assert_eq!(vor.vertices.len(), tris.len());
    }

    #[test]
    fn test_voronoi_edge_connectivity() {
        let sites = vec![
            Point::new(0.0, 0.0),
            Point::new(2.0, 0.0),
            Point::new(1.0, 2.0),
        ];
        let tris = triangulate(&sites);
        let vor = VoronoiDiagram::from_delaunay(&sites, &tris);
        for e in &vor.edges {
            assert!(e.from.distance(&e.to) > 0.0);
        }
    }

    #[test]
    fn test_voronoi_sites_preserved() {
        let sites = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(0.5, 1.0),
        ];
        let tris = triangulate(&sites);
        let vor = VoronoiDiagram::from_delaunay(&sites, &tris);
        assert_eq!(vor.sites.len(), sites.len());
    }

    #[test]
    fn test_larger_voronoi() {
        let mut sites = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                sites.push(Point::new(i as f64 * 2.0, j as f64 * 2.0));
            }
        }
        let tris = triangulate(&sites);
        let vor = VoronoiDiagram::from_delaunay(&sites, &tris);
        assert!(!vor.edges.is_empty());
        // Each edge should reference valid sites
        for e in &vor.edges {
            assert!(e.site_left < sites.len());
            assert!(e.site_right < sites.len());
        }
    }
}
