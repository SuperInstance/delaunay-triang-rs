//! Tutorial: Delaunay triangulation and Voronoi diagrams
//!
//! Shows Bowyer-Watson triangulation, edge flip optimization, and Voronoi extraction.

use delaunay_triang_rs::triangle::{Point, Triangle};
use delaunay_triang_rs::bowyer_watson::triangulate;
use delaunay_triang_rs::voronoi::VoronoiDiagram;

fn main() {
    println!("=== Delaunay Triangulation Tutorial ===\n");

    // Part 1: Basic triangulation
    println!("Part 1: Triangulating agent positions");
    let agents = vec![
        Point::new(0.0, 0.0),   // agent A
        Point::new(1.0, 0.0),   // agent B
        Point::new(0.5, 1.0),   // agent C
        Point::new(2.0, 0.5),   // agent D
        Point::new(1.5, 1.5),   // agent E
    ];
    
    let triangles = triangulate(&agents);
    println!("  {} agents → {} triangles", agents.len(), triangles.len());
    for (i, tri) in triangles.iter().enumerate() {
        println!("    Triangle {}: vertices {:?}", i, tri.vertices());
    }
    println!();

    // Part 2: Voronoi diagram from Delaunay
    println!("Part 2: Voronoi diagram (nearest-neighbor regions)");
    let voronoi = VoronoiDiagram::from_delaunay(&agents, &triangles);
    println!("  {} Voronoi edges", voronoi.edges.len());
    
    for i in 0..agents.len() {
        let region = voronoi.edges_for_site(i);
        println!("    Agent {} region: {} edges", i, region.len());
    }
    println!();

    // Part 3: Circumcircle test
    println!("Part 3: Point-in-circumcircle test (Delaunay criterion)");
    let a = Point::new(0.0, 0.0);
    let b = Point::new(1.0, 0.0);
    let c = Point::new(0.5, 0.866); // equilateral triangle
    
    let inside = Point::new(0.5, 0.3);
    let outside = Point::new(5.0, 5.0);
    
    println!("  Point ({}, {}) inside circumcircle: {}", inside.x, inside.y, inside.in_circumcircle(&a, &b, &c));
    println!("  Point ({}, {}) inside circumcircle: {}", outside.x, outside.y, outside.in_circumcircle(&a, &b, &c));
    println!();

    // Part 4: Distance calculations
    println!("Part 4: Distance calculations");
    let p1 = Point::new(0.0, 0.0);
    let p2 = Point::new(3.0, 4.0);
    println!("  Distance ({},{}) → ({},{}): {:.2}", p1.x, p1.y, p2.x, p2.y, p1.distance(&p2));
}
