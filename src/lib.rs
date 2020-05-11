//! # Rusty ABOS
//! 
//! Crate to generate a surface with the ABOS algorithm proposed in Art of Surface Interpolation by Mgr. Miroslav Dressler, Ph.D.
//!
//! http://m.dressler.sweb.cz/ABOS.htm
//!
//! Interpolates z values across a grid given a set of xyz points
//!
//! # Important Functions and Structs
//! Initalize `ABOSInputs` with points, set algorithm parameters, 
//! Use `abos_run_autogrid` to autoset grid parameters as specified by ABOS algorithm
//! Use `abos_run_grid_file` to use custom grid
//! Use `ABOSOutputs` to get desired output

extern crate approx;
extern crate nalgebra as na;

mod abos_constructor;
pub mod abos_structs;
mod io_system;

use crate::abos_constructor::{new_abos_autogrid, new_abos_grid_file};
use crate::abos_structs::{ABOSOutputs, ABOSInputs, ABOSImmutable,  ABOSMutable, INFINITY};
use crate::io_system::export_p_matrix;
use std::cmp;

/// Function to run abos, with algorithm recommended grid cell size, with extents by input points extents
pub fn abos_run_autogrid(abos_inputs: &ABOSInputs) -> ABOSOutputs{
    //1 Filtering points XYZ, specification of the grid, computation of the matrices NB and K, Z→DZ, 0→DP
    let (abos_immutable, mut abos_mutable) = new_abos_autogrid(&abos_inputs);
    abos_run(&mut abos_mutable,  &abos_immutable);
    ABOSOutputs::new(&abos_inputs, &abos_mutable, &abos_immutable)
}

/// Function to run abos, with grid file set cell size and extents
pub fn abos_run_grid_file(abos_inputs: &ABOSInputs, grid_file: String) -> ABOSOutputs{
    //1 Filtering points XYZ, specification of the grid, computation of the matrices NB and K, Z→DZ, 0→DP
    let (abos_immutable, mut abos_mutable) = new_abos_grid_file(&abos_inputs, grid_file);
    abos_run(&mut abos_mutable,  &abos_immutable);
    ABOSOutputs::new(&abos_inputs, &abos_mutable, &abos_immutable)
}

///Performs ABOS iteration cycle
/// 1. Filtering points XYZ, specification of the grid, computation of the matrices NB and K, Z→DZ, 0→DP
/// 2. Per partes constant interpolation of values DZ into the matrix P
/// 3. Tensioning, Linear Tensioning, and Smoothing of the matrix P
/// 4. P+DP→P
/// 5. Z - f(X Y) →DZi
/// 6. If the maximal difference max { DZ, ,i=1,..., n } does not exceed defined precision, the algorithm is finished
/// 7. P→DP, continue from step 2 again (= start the next iteration cycle)
pub fn abos_run(abos_mutable: &mut ABOSMutable, abos_immutable: &ABOSImmutable){
   // //Calculates k_u_v and kmax
   println!(
    "x cells {} y cells {} points {}",
    abos_immutable.i1,
    abos_immutable.j1,
    abos_immutable.xyz_points.len() / 3
    );
    //println!("--------//1 Initialization");
    //output_all_matrixes(&&abos_mutable, &abos_immutable);
    calculate_dz(abos_mutable, abos_immutable);

    let mut n = 1000;
    while n > 0 {
        n -= 1;
        //2 Per partes constant interpolation of values DZ into the matrix P
        per_parts_constant_interpolation(abos_mutable, abos_immutable);
        //export_p_matrix(&&abos_mutable, &abos_immutable, "afterPpli");
        //println!("--------//2 Per partes constant interpolation of values DZ into the matrix P");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);

        //3 Tensioning and smoothing of the matrix P
        tension_loop(abos_mutable, abos_immutable);

        //export_p_matrix(&&abos_mutable, &abos_immutable, "afterTension");
        //println!("--------//3 tensioning");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);

        linear_tension_loop(abos_mutable, abos_immutable);
        //export_p_matrix(&&abos_mutable, &abos_immutable, "afterLinearTension");

        //println!("--------//3 linear tensioning");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);

        smoothing_loop(abos_mutable, abos_immutable);
        //export_p_matrix(&&abos_mutable, &abos_immutable, "afterSmoothing");
        //println!("{:.1}", abos_mutable.p);
        //println!("--------//3 smoothing");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);

        //4 p + dp -> P
        abos_mutable.p += &abos_mutable.dp;
        //println!("--------//4 p + dp -> P");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);
        //back to 4

        //5 Dz - p => dz (Z - f (Xi, Yi) →DZi )
        calculate_dz(abos_mutable, abos_immutable);
        //println!("--------//5 Dz - p => dz (Z - f (Xi, Yi) →DZi )");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);
        // abos_mutable.dz -= abos_mutable.p;

        //6 calculate maximal difference
        //export_p_matrix(&&abos_mutable, &abos_immutable, "beforeDZCheck");
        let max_difference = abos_mutable.dz.max();
        println!("max_difference {}", max_difference);
        if max_difference.abs() < 0.000000001 {
            export_p_matrix(&&abos_mutable, &abos_immutable, "exportConverged");
            break;
        }
        //println!("--------//6 calculate maximal difference {}", max_difference);
        //output_all_matrixes(&&abos_mutable, &abos_immutable);

        //7 p -> Dp
        abos_mutable.dp.copy_from(&abos_mutable.p);
        //println!("--------//7 p -> Dp");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);
    }

    //output_all_matrixes(&&abos_mutable, &abos_immutable);
    //Initialization will swap so greater number of rows than cols
    if abos_immutable.xy_swaped {
        abos_mutable.p.swap_columns(0, 1);
    }
}

/// Calculate difference between predicted, and xyz point value in cell using bilinear interpolation.
fn calculate_dz(abos_mutable: &mut ABOSMutable, abos_immutable: &ABOSImmutable) {
    
    let mut dz_new = abos_immutable.z.clone();

    for (index, point) in abos_immutable.xyz_points.row_iter().enumerate() {
        let (x, y, _z) = (point[0], point[1], point[2]);
        let x_index_percentage = (x - abos_immutable.x1) / (abos_immutable.x2 - abos_immutable.x1); //contextualizes the index percentage
        let y_index_percentage = (y - abos_immutable.y1) / (abos_immutable.y2 - abos_immutable.y1);

        let x_index: f64 = (abos_immutable.i1 - 1) as f64 * x_index_percentage;
        let x_index_down: usize = x_index.floor() as usize;
        let mut x_index_up: usize = x_index.ceil() as usize;
        let value_index_percentage_x: f64 = x_index - x_index_down as f64;

        let y_index: f64 = (abos_immutable.j1 - 1) as f64 * y_index_percentage;
        let y_index_down: usize = y_index.floor() as usize;
        let mut y_index_up: usize = y_index_down + 1;
        let value_index_percentage_y: f64 = y_index - y_index_down as f64;

        if y_index_up == (abos_immutable.j1 as usize) {
            y_index_up -= 1;
        }
        if x_index_up == (abos_immutable.x1 as usize) {
            x_index_up -= 1;
        }
        let p = &abos_mutable.p;
        let mut corners = (0.0, 0.0, 0.0, 0.0);
        unsafe {
            corners.0 = *p.get_unchecked((x_index_down, y_index_down));
            corners.1 = *p.get_unchecked((x_index_down, y_index_up));
            corners.2 = *p.get_unchecked((x_index_up, y_index_down));
            corners.3 = *p.get_unchecked((x_index_up, y_index_up));
        }
        let x_down_val = ((1.0 - value_index_percentage_x) * corners.0)
            + ((value_index_percentage_x) * corners.1);
        let x_up_val = ((1.0 - value_index_percentage_x) * corners.2)
            + ((value_index_percentage_x) * corners.3);
        let y_val = ((1.0 - value_index_percentage_y) * x_down_val)
            + ((value_index_percentage_y) * x_up_val);

        unsafe {
            *dz_new.get_unchecked_mut(index) -= y_val;
        }
    }
    abos_mutable.dz.copy_from(&dz_new);
}

/// Calculating distance to reference cells in linear_tension_loop()
/// - See Art of Surface Interpolation 2.2.6
fn get_scaled_u_v(u: f64, v: f64, n: f64) -> (f64, f64) {
    let (mut u_mod, mut v_mod) = (u, v);
    let uv_magnitude = (u * u + v * v).sqrt();
    if uv_magnitude > n {
        let c = n / uv_magnitude;
        u_mod = c * u;
        v_mod = c * v;
    };
    (u_mod, v_mod)
}

/// Makes a weighted average of 4 corner cells, top right, top left, bot right, bot left
/// - See Art of Surface Interpolation 2.2.5
//TODO shift from corner check, to mid perimeter check
fn tension_cell(ii: i32, jj: i32, k_i_j_mod: i32, abos_mutable: &mut ABOSMutable) -> f64 {
    //we need to get Pi , j=Pik , jPi , jkPi−k , jPi , j−k
    let p = &abos_mutable.p;

    let (min_i, max_i) = get_valid_dim_bounds(ii as usize, k_i_j_mod as usize, 0, abos_mutable.p.nrows());
    let (min_j, max_j) = get_valid_dim_bounds(jj as usize, k_i_j_mod as usize, 0, abos_mutable.p.ncols());

    let mut p1: f64 = 0.0;
    unsafe {
        p1 += (*p.get_unchecked((min_i , min_j)))
            + (*p.get_unchecked((max_i,  max_j)))
            + (*p.get_unchecked((min_i,  max_j)))
            + (*p.get_unchecked((max_i,  min_j)));
        p1 /= 4.0;
    }
    p1
}

/// Sets Z value of cell to nearest known value, equivalent to voronoi
/// - See Art of Surface Interpolation 2.2.4
pub fn per_parts_constant_interpolation(
    abos_mutable: &mut ABOSMutable,
    abos_immutable: &ABOSImmutable,
) {
    for (ii, row) in abos_immutable.nb.row_iter().enumerate() {
        for (jj, col) in row.iter().enumerate() {
            let dz_closest: f64 = *abos_mutable.dz.get(*col).unwrap();
            unsafe {
                let p_position = abos_mutable.p.get_unchecked_mut((ii, jj));
                *p_position = dz_closest;
            }
        }
    }
}

/// Makes a weighted average of 4 corner cells, top right, top left, bot right, bot left
/// - See Art of Surface Interpolation 2.2.5
pub fn tension_loop(abos_mutable: &mut ABOSMutable, abos_immutable: &ABOSImmutable) {
    let n = cmp::max(4, abos_immutable.k_max / 2 + 2); //TODO SWITCH BACK
    for n_countdown in (0..n).rev() {
        for (ii, row) in abos_immutable.k_u_v.row_iter().enumerate() {
            for (jj, col) in row.iter().enumerate() {
                let mut k_i_j_mod = col.0; //If  k  is greater than the decreasing loop variable N, then k = N.
                if k_i_j_mod > (n_countdown + 1) {
                    k_i_j_mod = (n_countdown + 1);
                }
                let new_p = tension_cell(ii as i32, jj as i32, k_i_j_mod as i32, abos_mutable);
                unsafe {
                    if new_p != -INFINITY {
                        *abos_mutable.p.get_unchecked_mut((ii, jj)) = new_p;
                    }
                }
            }
        }
    }
}

/// Calculates modified weighted average value for cell
/// - See Art of Surface Interpolation 2.2.6
/// - Has considerations and weighting for corner and edge cell not specified by paper
fn linear_tension_cell(
    q: f64,
    ii: usize,
    jj: usize,
    u: usize,
    v: usize,
    abos_immutable: &ABOSImmutable,
    abos_mutable: &mut ABOSMutable,
) -> f64 {
    let p = &abos_mutable.p;
    let (top_l, top_l_b) = if (ii + u) >= p.nrows() as usize || (jj + v) >= p.ncols() {
        (0.0, 0.0)
    } else {
        unsafe { (*p.get_unchecked((ii + u, jj + v)), 1.0) }
    };

    let (bot_r, bot_r_b) = if (ii as i32 - u as i32) < 0 || (jj as i32 - v as i32) < 0 {
        (0.0, 0.0)
    } else {
        unsafe { (*p.get_unchecked((ii - u, jj - v)), 1.0) }
    };

    let (top_r, top_r_b) = if (ii as i32 - u as i32) < 0 || jj + v >= p.ncols() {
        (0.0, 0.0)
    } else {
        unsafe { (*p.get_unchecked((ii - u, jj + v)), 1.0) }
    };

    let (bot_l, bot_l_b) = if ii + u >= p.nrows() as usize || (jj as i32 - v as i32) < 0 {
        (0.0, 0.0)
    } else {
        unsafe { (*p.get_unchecked((ii + u, jj - v)), 1.0) }
    };
    let bottom_divider = (q * top_l_b)
        + (q * bot_r_b)
        + top_r_b * abos_immutable.r as f64
        + bot_l_b * abos_immutable.r as f64;

    if bottom_divider == 0.0 {
        -INFINITY
    } else {
        (q * top_l + q * bot_r + top_r * abos_immutable.r as f64 + bot_l * abos_immutable.r as f64)
            / (bottom_divider)
    }
}

/// Modifies the matrix P according to the formula for weighted average
/// - Has degrees of linear tension specified by user inputs
/// - See Art of Surface Interpolation 2.2.6
pub fn linear_tension_loop(abos_mutable: &mut ABOSMutable, abos_immutable: &ABOSImmutable) {
    let n = cmp::min(4, abos_immutable.k_max / 2 + 2); //TODO SWITCH BACK
    for n_countdown in (1..n + 1).rev() {
        for (ii, row) in abos_immutable.k_u_v.row_iter().enumerate() {
            for (jj, col) in row.iter().enumerate() {
                let k_i_j_mod = cmp::min(col.0, n_countdown);
                let q = get_q_linear_tension(abos_immutable, k_i_j_mod);
                let (u_mod, v_mod) = get_scaled_u_v(col.1 as f64, col.2 as f64, n as f64);
                let new_p = linear_tension_cell(
                    q,
                    ii,
                    jj,
                    u_mod as usize,
                    v_mod as usize,
                    abos_immutable,
                    abos_mutable,
                );
                if new_p != -INFINITY {
                    unsafe {
                        *abos_mutable.p.get_unchecked_mut((ii, jj)) = new_p;
                    }
                }
            }
        }
    }
}

/// Adjustemnts to linear tension specified in
/// - See Art of Surface Interpolation 3.4.2 Degrees of linear tensioning.
fn get_q_linear_tension(abos_immutable: &ABOSImmutable, k_i_j: usize) -> f64 {
    //println!("self.l {}, self.k_max {}", self.l, self.k_max);
    match abos_immutable.degree {
        3 => 1.0,
        2 => abos_immutable.l * (abos_immutable.k_max - k_i_j) as f64,
        _ => abos_immutable.l * ((abos_immutable.k_max - k_i_j) as f64).powi(2),
    }
}

/// Function to manage bound checking with edge and corner cells
/// - tt : midpoint, checking window from
/// - dt : amount to look above below
/// - min_i : minimum acceptable value, inclusive
/// - max_i  : maximum acceptable value exclusive
/// - returns (min_bound, max_bound)  both inclusive indexes
/// - where min_bound = tt-dt || min_t
/// - where max_bound = tt+dt || max_t - 1
pub fn get_valid_dim_bounds(tt: usize, dt: usize, min_t: usize, max_t: usize) -> (usize, usize) {
    let lower_bound = if tt < dt || tt - dt < min_t {
            min_t
        } else {
            tt - dt
        };
    let upper_bound = if tt + dt >= max_t { 
            max_t - 1
        } else { 
            tt + dt 
        };
    (lower_bound, upper_bound)
}

/// Setting smoothing parameter for each cell
/// Art of Surface Interpolation 2.2.7 Smoothing.
pub fn set_t_smooth(abos_mutable: &mut ABOSMutable) {
    //Set initial tt matrix
    for (ii, row) in abos_mutable.p.row_iter().enumerate() {
        let (kk_min, kk_max) = get_valid_dim_bounds(ii, 2, 0, abos_mutable.p.nrows());
        for (jj, _col) in row.iter().enumerate() {
            let (ll_min, ll_max) = get_valid_dim_bounds(jj, 2, 0, abos_mutable.p.ncols());
            let num_cells = ((1 + kk_max - kk_min) * (1 + ll_max - ll_min)) as f64;
            //println!("ii {}, jj {},  kk_min {}, kk_max {}, ll_min {}, ll_max{}, nrows {}, ncols {}", ii, jj, kk_min, kk_max, ll_min, ll_max, abos_mutable.p.nrows(), abos_mutable.p.ncols());
            //.slice(start, shape)
            let mut pij_resid_sum = abos_mutable
                .p
                .slice((kk_min, ll_min), (kk_max - kk_min + 1, ll_max - ll_min + 1)) //(kk_max - kk_min + 1, ll_max - ll_min + 1)
                .sum();
            unsafe {
                pij_resid_sum -= num_cells * (*abos_mutable.p.get_unchecked((ii, jj)));
                *abos_mutable.t_smooth.get_unchecked_mut((ii, jj)) = pij_resid_sum * pij_resid_sum;
            }
        }
    }

    //scale t_smooth matrix to between 0 to 1.0, variation on published method
    let min_t = abos_mutable.t_smooth.min();
    let dt = abos_mutable.t_smooth.max() - min_t;
    if dt > 0.0  {
        abos_mutable.t_smooth.apply(|x| (x - min_t) * 1.0 / dt);
    }
}

/// Setting smoothing parameter for each cell
/// Has degrees of linear tension specified by user inputs
/// See Art of Surface Interpolation 2.2.7 Smoothing
pub fn smoothing_loop(abos_mutable: &mut ABOSMutable, abos_immutable: &ABOSImmutable) {
    let n = cmp::max(4, abos_immutable.k_max * (abos_immutable.k_max / 16));

    //t , are weights, which are zero before the first smoothing and afterwards they are computed according to the formula
    // abos_mutable.t_smooth.fill(0.0);
    set_t_smooth(abos_mutable);
    for n_countdown in (1..n + 1).rev() {
        for (ii, row) in abos_immutable.k_u_v.row_iter().enumerate() {
            for (jj, col) in row.iter().enumerate() {
                let k_i_j = col.0;
                if k_i_j >= n_countdown {
                    let (kk_min, kk_max) =
                        get_valid_dim_bounds(ii, 1, 0, abos_mutable.p.nrows());
                    let (ll_min, ll_max) =
                        get_valid_dim_bounds(jj, 1, 0, abos_mutable.p.ncols());

                    let num_cells = ((kk_max - kk_min + 1) * (ll_max - ll_min + 1)) as f64; //TODO THIS IS DECREMENTING THE TOTAL VALUE OF THE GRID

                    let mut p_i_j_val: f64 = -INFINITY;
                    let mut t_i_j_val: f64 = -INFINITY;
                    unsafe {
                        p_i_j_val = *abos_mutable.p.get_unchecked((ii, jj));
                        t_i_j_val = *abos_mutable.t_smooth.get_unchecked((ii, jj));
                    }
                    let sum_p_k_l = abos_mutable
                        .p
                        .slice((kk_min, ll_min), (kk_max - kk_min + 1, ll_max - ll_min + 1))
                        .sum();
                    let right_modified = p_i_j_val * (abos_immutable.q_smooth * t_i_j_val - 1.0);

                    let new_p_i_j = (sum_p_k_l + right_modified)
                        / (num_cells - 1.0 + t_i_j_val * abos_immutable.q_smooth);

                    unsafe {
                        *abos_mutable.p.get_unchecked_mut((ii, jj)) = new_p_i_j;
                    }
                }
            }
        }
    }
}

fn output_all_matrixes(abos_mutable: &ABOSMutable, abos_immutable: &ABOSImmutable) {
    println!(
        "NB {:.1} Z{:.1} DZ{:.1} DP{:.1} P{:.1}",
        abos_immutable.nb, abos_immutable.z, abos_mutable.dz, abos_mutable.dp, abos_mutable.p
    );
    //println!("{}", abos_mutable.p.len())
}