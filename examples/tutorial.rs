//! # delaunay-triang-rs Tutorial
//!
//! A progressive walkthrough of Delaunay triangulation, Voronoi diagrams,
//! edge-flip optimization, and quad-edge data structures in pure Rust.
//!
//! ## Lessons
//!
//! 1. **Points & Triangles** — primitives, distances, circumcircles, area
//! 2. **Bowyer-Watson Triangulation** — compute Delaunay triangulation from points
//! 3. **Voronoi Diagrams** — the dual of Delaunay, querying edges per site
//! 4. **Edge-Flip Algorithm** — convert any triangulation into Delaunay
//! 5. **Quad-Edge Subdivision** — Guibas-Stolfi topological data structure
//! 6. **Putting It All Together** — mesh a grid, compute Voronoi, verify quality

fn main() {
    println!("=== delaunay-triang-rs Tutorial ===\n");

    lesson_1_points_and_triangles();
    lesson_2_bowyer_watson();
    lesson_3_voronoi_diagrams();
    lesson_4_edge_flip();
    lesson_5_quad_edge();
    lesson_6_end_to_end();
}

// ─── Lesson 1: Points & Triangles ────────────────────────────────────────────

fn lesson_1_points_and_triangles() {
    println!("Lesson 1: Points & Triangles");
    println!("-----------------------------\n");

    use delaunay_triang_rs::{Point, Triangle};

    // --- Creating points ---
    let a = Point::new(0.0, 0.0);
    let b = Point::new(3.0, 0.0);
    let c = Point::new(0.0, 4.0);
    println!("Points: a={}, b={}, c={}", a, b, c);

    // --- Distance computation ---
    let dist_ab = a.distance(&b);
    let dist_sq_ab = a.dist_sq(&b);
    println!(
        "Distance a→b = {:.4}  (squared = {:.4})",
        dist_ab, dist_sq_ab
    );
    assert!((dist_ab - 3.0).abs() < 1e-10);

    // --- Triangles from vertex indices ---
    let tri = Triangle::new(0, 1, 2);
    println!("Triangle vertices: {:?}", tri.vertices());
    assert!(tri.contains_vertex(1));
    assert!(!tri.contains_vertex(5));

    // --- Circumcircle ---
    let points = [a, b, c];
    let (center, r_sq) = tri.circumcircle(&points);
    println!("Circumcircle center = {}, radius² = {:.4}", center, r_sq);

    // All vertices are equidistant from the circumcenter
    assert!((center.dist_sq(&points[0]) - r_sq).abs() < 1e-6);
    assert!((center.dist_sq(&points[1]) - r_sq).abs() < 1e-6);
    assert!((center.dist_sq(&points[2]) - r_sq).abs() < 1e-6);

    // --- Area ---
    let area = tri.area(&points);
    let signed = tri.signed_area(&points);
    println!("Area = {:.4}, Signed area = {:.4}", area, signed);
    assert!((area - 6.0).abs() < 1e-10); // 3-4-5 triangle → area = 6

    // --- In-circumcircle test ---
    let interior = Point::new(0.5, 0.5);
    let exterior = Point::new(10.0, 10.0);
    assert!(interior.in_circumcircle(&a, &b, &c));
    assert!(!exterior.in_circumcircle(&a, &b, &c));
    println!("Interior point in circumcircle: true");
    println!("Exterior point in circumcircle: false");

    // --- Edge sharing ---
    let other = Triangle::new(1, 2, 3);
    let isolated = Triangle::new(4, 5, 6);
    assert!(tri.shares_edge(&other));
    assert!(!tri.shares_edge(&isolated));
    println!("Triangle shares edge with [1,2,3]: true");
    println!("Triangle shares edge with [4,5,6]: false");

    println!();
}

// ─── Lesson 2: Bowyer-Watson Triangulation ────────────────────────────────────

fn lesson_2_bowyer_watson() {
    println!("Lesson 2: Bowyer-Watson Triangulation");
    println!("--------------------------------------\n");

    use delaunay_triang_rs::bowyer_watson::triangulate;
    use delaunay_triang_rs::{Point, Triangle};

    // --- Triangulate a simple square ---
    let square = vec![
        Point::new(0.0, 0.0),
        Point::new(1.0, 0.0),
        Point::new(1.0, 1.0),
        Point::new(0.0, 1.0),
    ];

    let tris: Vec<Triangle> = triangulate(&square);
    println!("Square (4 points) → {} triangles", tris.len());
    assert_eq!(tris.len(), 2);

    for (i, t) in tris.iter().enumerate() {
        println!(
            "  Triangle {}: vertices = {:?}, area = {:.4}",
            i,
            t.vertices(),
            t.area(&square)
        );
    }

    // Total area should equal 1.0 (the unit square)
    let total_area: f64 = tris.iter().map(|t| t.area(&square)).sum();
    assert!((total_area - 1.0).abs() < 1e-6);
    println!("  Total area = {:.6}", total_area);

    // --- Fewer than 3 points returns empty ---
    let few = vec![Point::new(0.0, 0.0), Point::new(1.0, 0.0)];
    assert!(triangulate(&few).is_empty());
    println!("  <3 points → 0 triangles (expected)");

    // --- Verify the Delaunay property on a larger set ---
    let mut pts = Vec::new();
    for i in 0..5 {
        for j in 0..5 {
            pts.push(Point::new(i as f64 * 2.0, j as f64 * 2.0));
        }
    }
    let tris = triangulate(&pts);
    println!("5×5 grid (25 points) → {} triangles", tris.len());

    // Delaunay invariant: no point lies inside any triangle's circumcircle
    for t in &tris {
        let (center, r_sq) = t.circumcircle(&pts);
        for (i, p) in pts.iter().enumerate() {
            if t.contains_vertex(i) {
                continue;
            }
            assert!(
                center.dist_sq(p) >= r_sq - 1e-6,
                "Delaunay property violated at point {}",
                i
            );
        }
    }
    println!("  ✓ Delaunay property verified for all triangles");

    println!();
}

// ─── Lesson 3: Voronoi Diagrams ──────────────────────────────────────────────

fn lesson_3_voronoi_diagrams() {
    println!("Lesson 3: Voronoi Diagrams");
    println!("---------------------------\n");

    use delaunay_triang_rs::bowyer_watson::triangulate;
    use delaunay_triang_rs::voronoi::{VoronoiDiagram, VoronoiEdge};
    use delaunay_triang_rs::Point;

    // --- Build Voronoi from Delaunay ---
    let sites = vec![
        Point::new(0.0, 0.0),
        Point::new(4.0, 0.0),
        Point::new(2.0, 3.0),
        Point::new(2.0, -1.0),
    ];

    let delaunay = triangulate(&sites);
    let voronoi = VoronoiDiagram::from_delaunay(&sites, &delaunay);

    println!("Sites: {}", sites.len());
    println!("Voronoi vertices: {}", voronoi.vertices.len());
    println!("Voronoi edges: {}", voronoi.edges.len());

    // --- Inspect edges ---
    for (i, edge) in voronoi.edges.iter().enumerate() {
        println!(
            "  Edge {}: ({:.2},{:.2}) → ({:.2},{:.2})  [sites {} | {}]",
            i, edge.from.x, edge.from.y, edge.to.x, edge.to.y, edge.site_left, edge.site_right
        );
    }

    // --- Query edges for a specific site ---
    let site_0_edges: Vec<&VoronoiEdge> = voronoi.edges_for_site(0);
    println!(
        "  Site 0 has {} neighboring Voronoi edges",
        site_0_edges.len()
    );

    // --- Larger grid ---
    let mut grid = Vec::new();
    for i in 0..4i32 {
        for j in 0..4i32 {
            grid.push(Point::new(i as f64, j as f64));
        }
    }
    let d = triangulate(&grid);
    let v = VoronoiDiagram::from_delaunay(&grid, &d);
    println!(
        "4×4 grid → {} Voronoi edges, {} vertices",
        v.edges.len(),
        v.vertices.len()
    );

    // All edge endpoints should reference valid sites
    for e in &v.edges {
        assert!(e.site_left < grid.len());
        assert!(e.site_right < grid.len());
    }
    println!("  ✓ All edge site references are valid");

    println!();
}

// ─── Lesson 4: Edge-Flip Algorithm ───────────────────────────────────────────

fn lesson_4_edge_flip() {
    println!("Lesson 4: Edge-Flip Algorithm");
    println!("------------------------------\n");

    use delaunay_triang_rs::edge_flip::edge_flip_delaunay;
    use delaunay_triang_rs::{Point, Triangle};

    // --- Start with a non-Delaunay triangulation and fix it ---
    let pts = vec![
        Point::new(0.0, 0.0),
        Point::new(2.0, 0.0),
        Point::new(1.0, 2.0),
        Point::new(1.0, 0.5), // interior point
    ];

    // Manually create a non-optimal triangulation
    let mut tris = vec![
        Triangle::new(0, 1, 3),
        Triangle::new(0, 2, 3),
        Triangle::new(1, 2, 3),
    ];

    println!("Before edge flip: {} triangles", tris.len());
    for (i, t) in tris.iter().enumerate() {
        println!(
            "  T{}: vertices {:?}, area = {:.4}",
            i,
            t.vertices(),
            t.area(&pts)
        );
    }

    // Apply edge flips to enforce the Delaunay property
    edge_flip_delaunay(&pts, &mut tris);

    println!("\nAfter edge flip: {} triangles", tris.len());
    for (i, t) in tris.iter().enumerate() {
        println!(
            "  T{}: vertices {:?}, area = {:.4}",
            i,
            t.vertices(),
            t.area(&pts)
        );
    }

    // Verify Delaunay property
    for t in &tris {
        let (center, r_sq) = t.circumcircle(&pts);
        for (i, p) in pts.iter().enumerate() {
            if t.contains_vertex(i) {
                continue;
            }
            assert!(
                center.dist_sq(p) >= r_sq - 1e-4,
                "Post-flip Delaunay violation"
            );
        }
    }
    println!("  ✓ Delaunay property satisfied after flipping");

    // --- Square: already Delaunay, no flips needed ---
    let sq_pts = vec![
        Point::new(0.0, 0.0),
        Point::new(1.0, 0.0),
        Point::new(1.0, 1.0),
        Point::new(0.0, 1.0),
    ];
    let mut sq_tris = vec![Triangle::new(0, 1, 2), Triangle::new(0, 2, 3)];
    let count_before = sq_tris.len();
    edge_flip_delaunay(&sq_pts, &mut sq_tris);
    assert_eq!(sq_tris.len(), count_before);
    println!("  Square: no flip needed (already Delaunay)");

    println!();
}

// ─── Lesson 5: Quad-Edge Data Structure ──────────────────────────────────────

fn lesson_5_quad_edge() {
    println!("Lesson 5: Quad-Edge Data Structure");
    println!("------------------------------------\n");

    use delaunay_triang_rs::quad_edge::QuadEdgeSubdivision;
    use delaunay_triang_rs::Point;

    // --- Build a triangle from three edges ---
    let mut qes = QuadEdgeSubdivision::new();
    let a = qes.add_point(Point::new(0.0, 0.0));
    let b = qes.add_point(Point::new(1.0, 0.0));
    let c = qes.add_point(Point::new(0.5, 0.866));

    println!("Added 3 points: indices {}, {}, {}", a, b, c);

    // Create edges forming a triangle
    let e_ab = qes.make_edge(a, b);
    let e_bc = qes.make_edge(b, c);
    let e_ca = qes.make_edge(c, a);

    // Splice edges to form a closed loop
    qes.splice(e_ab.sym(), e_bc);
    qes.splice(e_bc.sym(), e_ca);
    qes.splice(e_ca.sym(), e_ab);

    println!("Created triangle with {} edges", qes.edge_count());
    assert_eq!(qes.edge_count(), 3);

    // --- Inspect edge topology ---
    assert_eq!(qes.origin(e_ab), Some(a));
    assert_eq!(qes.dest(e_ab), Some(b));
    println!("Edge a→b: origin={:?}, dest={:?}", qes.origin(e_ab), qes.dest(e_ab));

    // Symmetric edge reverses direction
    let e_ba = e_ab.sym();
    assert_eq!(qes.origin(e_ba), Some(b));
    assert_eq!(qes.dest(e_ba), Some(a));
    println!("Edge b→a (sym): origin={:?}, dest={:?}", qes.origin(e_ba), qes.dest(e_ba));

    // --- Edge rotations ---
    // Rotating twice gives the symmetric edge
    assert_eq!(e_ab.rot().rot(), e_ab.sym());
    println!("rot(rot(e)) == sym(e): verified");

    // --- Delete an edge ---
    let mut qes2 = QuadEdgeSubdivision::new();
    let p0 = qes2.add_point(Point::new(0.0, 0.0));
    let p1 = qes2.add_point(Point::new(1.0, 0.0));
    let e = qes2.make_edge(p0, p1);
    assert_eq!(qes2.edge_count(), 1);
    qes2.delete_edge(e);
    println!("Edge deleted (topological disconnect)");

    // --- Build a larger subdivision: 4-cycle ---
    let mut ring = QuadEdgeSubdivision::new();
    let pts: Vec<usize> = (0..4)
        .map(|i| ring.add_point(Point::new((i as f64).cos(), (i as f64).sin())))
        .collect();

    for i in 0..4 {
        ring.make_edge(pts[i], pts[(i + 1) % 4]);
    }
    println!("4-cycle subdivision: {} edges", ring.edge_count());
    assert_eq!(ring.edge_count(), 4);

    println!();
}

// ─── Lesson 6: End-to-End Pipeline ───────────────────────────────────────────

fn lesson_6_end_to_end() {
    println!("Lesson 6: End-to-End Pipeline");
    println!("------------------------------\n");

    use delaunay_triang_rs::bowyer_watson::triangulate;
    use delaunay_triang_rs::voronoi::VoronoiDiagram;
    use delaunay_triang_rs::{Point, Triangle};

    // --- Generate random-ish point cloud ---
    let pts: Vec<Point> = (0..20)
        .map(|i| {
            let angle = i as f64 * 0.314159; // ~golden-angle-ish
            let r = 2.0 + (i as f64 * 0.1).sin() * 1.5;
            Point::new(r * angle.cos(), r * angle.sin())
        })
        .collect();
    println!("Generated {} points in a spiral pattern", pts.len());

    // --- Delaunay triangulation ---
    let tris: Vec<Triangle> = triangulate(&pts);
    println!("Delaunay: {} triangles", tris.len());

    // Verify total area is positive and all triangles are valid
    let total_area: f64 = tris.iter().map(|t| t.area(&pts)).sum();
    println!("Total mesh area = {:.4}", total_area);
    assert!(total_area > 0.0);

    for t in &tris {
        assert!(t.area(&pts) > 1e-12, "Degenerate triangle found");
    }
    println!("  ✓ All triangles have positive area");

    // --- Voronoi diagram ---
    let voronoi = VoronoiDiagram::from_delaunay(&pts, &tris);
    println!(
        "Voronoi: {} edges, {} vertices",
        voronoi.edges.len(),
        voronoi.vertices.len()
    );

    // Each interior site should have at least one edge
    let mut sites_with_edges = 0;
    for i in 0..pts.len() {
        if !voronoi.edges_for_site(i).is_empty() {
            sites_with_edges += 1;
        }
    }
    println!(
        "  {}/{} sites have Voronoi edges",
        sites_with_edges, pts.len()
    );

    // --- Verify Delaunay invariant ---
    for t in &tris {
        let (center, r_sq) = t.circumcircle(&pts);
        for (i, p) in pts.iter().enumerate() {
            if t.contains_vertex(i) {
                continue;
            }
            assert!(
                center.dist_sq(p) >= r_sq - 1e-6,
                "Invariant violated at point {}",
                i
            );
        }
    }
    println!("  ✓ Delaunay invariant holds across entire mesh");

    // --- Compute per-triangle quality (aspect ratio proxy) ---
    let mut min_quality = f64::MAX;
    for t in &tris {
        let a = t.area(&pts);
        let (center, r_sq) = t.circumcircle(&pts);
        // Quality metric: 4√3 · area / perimeter² (1.0 = equilateral)
        let quality = if r_sq > 1e-12 {
            4.0 * a / (r_sq * 3.0 * 1.7320508)
        } else {
            0.0
        };
        min_quality = min_quality.min(quality);
    }
    println!("  Worst triangle quality: {:.4}", min_quality);

    println!("\n=== Tutorial Complete! ===");
}
