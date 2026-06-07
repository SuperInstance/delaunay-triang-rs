//! # delaunay-triang-rs
//!
//! Delaunay triangulation and Voronoi diagram algorithms in pure Rust.
//!
//! # Algorithms
//!
//! - **Bowyer-Watson** — Incremental point insertion for Delaunay triangulation
//! - **Voronoi** — Dual graph computation from Delaunay triangulation
//! - **Edge Flip** — Convert any triangulation to Delaunay
//! - **Quad-Edge** — Guibas-Stolfi data structure

pub mod triangle;
pub mod bowyer_watson;
pub mod voronoi;
pub mod edge_flip;
pub mod quad_edge;

pub use triangle::{Point, Triangle};
