
extern crate approx;
extern crate nalgebra as na;

const INFINITY: f64 = 1.0f64 / 0.0f64;
use kdtree::KdTree;
use kdtree::distance::squared_euclidean;
use na::{DVector, MatrixMN, Dynamic, RealField, U3, U1};
extern crate rand;
use rand::prelude::*;
use crate::{initialize_dmatrix, ABOSGrid, tension_cell};
    

pub fn get_min_chebyshev_distance_n2(xyz_points: &MatrixMN<f64, Dynamic, U3>) -> f64 {
    let mut min_chebyshev_distance: f64 = INFINITY;

    //step 1: make a 2d search tree to accelerate the finding of points
    let dimensions = 2;
    let mut kdtree = KdTree::new(dimensions);
    // let mut points:Vec<([f64; 2], usize)> = vec![];
    for (i, row) in xyz_points.row_iter().enumerate() {
        kdtree.add([row[0], row[1]], i).unwrap();
    }

    let xyz_iter1 = xyz_points.row_iter();
    for (_, row) in xyz_iter1.enumerate() {
        let xyz_iter2 = xyz_points.row_iter();
        for (_, row2) in xyz_iter2.enumerate() {
            if row != row2 {
                let distances: MatrixMN<f64, U1, U3> = row.clone_owned() - row2;
                let distances = distances.abs();
                let max_xy_distance = if distances[0] > distances[1] { distances[0] } else { distances[1] };

                if max_xy_distance < min_chebyshev_distance {
                    min_chebyshev_distance = max_xy_distance;
                }
            }
        }
    }
    return min_chebyshev_distance;
}

fn get_min_chebyshev_distance_kd(xyz_points: &MatrixMN<f64, Dynamic, U3>) -> f64 {
    let mut min_chebyshev_distance: f64 = INFINITY;

    //step 1: make a 2d search tree to accelerate the finding of points
    let dimensions = 2;
    let mut kdtree = KdTree::new(dimensions);
    // let mut points:Vec<([f64; 2], usize)> = vec![];
    for (i, row) in xyz_points.row_iter().enumerate() {
        kdtree.add([row[0], row[1]], i).unwrap();
    }
    //TODO impliment n nearest because kdtree includesitself therefore nearest point is self
    for row in xyz_points.row_iter() {
        let point = [row[0], row[1]];
        let kd_search_result = kdtree.nearest(&point, 2, &squared_euclidean).unwrap();
        let closest_indx = *kd_search_result[1].1;
        let distances: MatrixMN<f64, U1, U3> = row.clone_owned() - xyz_points.row(closest_indx);
        let distances = distances.abs();
        let max_xy_distance = if distances[0] > distances[1] { distances[0] } else { distances[1] };

        if max_xy_distance < min_chebyshev_distance {
            min_chebyshev_distance = max_xy_distance;
        }
    }
    return min_chebyshev_distance;
}

pub fn indexes_to_position_direct<N: RealField>(indxs: DVector<N>, d_mins: DVector<N>, d_sizes: DVector<N>) -> DVector<N> {
    let position = d_mins + d_sizes * indxs;
    position
}

pub fn indexes_to_position(indxs: DVector<u32>, d_mins: DVector<f64>, d_sizes: DVector<f64>) -> DVector<f64> {
    let indxs_na: DVector<f64> = na::convert(indxs);
    indexes_to_position_direct(indxs_na, d_mins, d_sizes)
}


pub fn test_min_chebyshev_dist() {
    //-------- Testing Veracity of kd  --------------
    for _ in 0..500 {
        let mut points: Vec<Vec<f64>> = vec![];
        let mut rng = rand::thread_rng();
        for _ in 0..500 {
            let new_point:Vec<f64> = vec!(rng.gen_range(200.0,300.0),rng.gen_range(100.0,200.0),rng.gen_range(0.0,10.0));
            points.push(new_point);
        }
        let xyz_points: MatrixMN<f64, Dynamic, U3> = initialize_dmatrix(points);
        assert_eq!(
            get_min_chebyshev_distance_n2(&xyz_points),
            get_min_chebyshev_distance_kd(&xyz_points)
        );
    }
}


//justpass it an ABOSGrid parameter... vs passing 5 parameters..., 
pub fn tension_loop_functional(n: usize, i1: usize, j1: usize, mutable_p: &mut MatrixMN<f64, Dynamic, Dynamic>, 
    k: &MatrixMN<(usize, usize, usize), Dynamic, Dynamic>) {
    for n_countdown in (1..n + 1).rev() {
        for (ii, row) in k.row_iter().enumerate() {
            for (jj, col) in row.iter().enumerate() {
                let k_i_j_mod = std::cmp::min(col.0, n_countdown);
                let new_p = tension_cell(ii as i32, jj as i32, k_i_j_mod as i32, mutable_p);
                unsafe {
                    *mutable_p.get_unchecked_mut((ii, jj)) = new_p;
                }
            }
        }
    }
}

#[test]
pub fn test_tension_loop() {
    let mut points: Vec<Vec<f64>> = vec![];
    let mut points2: Vec<Vec<f64>> = vec![];
    for ii in 0..3 {
        let iif = 1.0 + ii as f64;
        let new_point:Vec<f64> = vec!(iif , iif*2.0, iif*3.0);
        let new_point2:Vec<f64> = vec!(iif , iif*2.0, iif*3.0);
        points.push(new_point);
        points2.push(new_point2);
    }
    
    let mut agrid1 = ABOSGrid::new(points, 30.0, 0);
    let mut agrid2 = ABOSGrid::new(points2, 30.0, 0);
    assert_eq!(agrid1.p, agrid2.p);
    assert_eq!(agrid1.k_u_v, agrid2.k_u_v);
    
    agrid1.per_parts_constant_interpolation();
    agrid2.per_parts_constant_interpolation();
    assert_eq!(agrid1.p, agrid2.p);
    assert_eq!(agrid1.k_u_v, agrid2.k_u_v);

    tension_loop_functional(agrid1.n, agrid1.i1 as usize, agrid1.j1 as usize, &mut agrid1.p, & agrid1.k_u_v);
    ABOSGrid::tension_loop(&mut agrid2);        
    assert_eq!(agrid1.p, agrid2.p);
    assert_eq!(agrid1.k_u_v, agrid2.k_u_v);
}