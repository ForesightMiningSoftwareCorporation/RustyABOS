///Crate to generate a surface with the ABOS algorithm proposed in Art of Surface Interpolation byM gr. Miroslav Dressler, Ph.D.
/// http://m.dressler.sweb.cz/ABOS.htm
///
extern crate approx;
extern crate nalgebra as na;
mod abos_constructor;
pub mod abos_structs;
mod io_system;
use crate::abos_constructor::new_abos;
use crate::abos_structs::{ABOSImmutable, ABOSInputs, ABOSMutable, INFINITY};
use na::{DVector, Dim, Dynamic, MatrixMN, U1, U3};
use std::cmp;

//iteration cycle
// 1. Filtering points XYZ, specification of the grid, computation of the matrices NB and
// K, Z→DZ, 0→DP
//test.init_distance_point_matrixes(); //prepare nb/k/kmax
// 2. Per partes constant interpolation of values DZ into the matrix P
// 3. Tensioning and smoothing of the matrix P
//----3 sub steps Tensioning, Linear Tensioning, Smoothing
// 4. P+DP→P
// 5. Z - f(X Y) →DZi
// 6. If the maximal difference max { DZ, ,i=1,..., n } does not exceed defined precision, the algorithm is finished
// 7. P→DP, continue from step 2 again (= start the next iteration cycle)

pub fn abos_run(abos_inputs: &ABOSInputs) {
    //1 Initialization
    let (abos_immutable, mut abos_mutable) = new_abos(&abos_inputs);
    // //Calculates k_u_v and kmax
    println!(
        "x cells {} y cells {}",
        abos_immutable.i1, abos_immutable.j1
    );
    println!("--------//1 Initialization");
    output_all_matrixes(&&abos_mutable, &abos_immutable);

    loop {
        //2 Per partes constant interpolation of values DZ into the matrix P
        per_parts_constant_interpolation(&mut abos_mutable, &abos_immutable);
        //println!("--------//2 Per partes constant interpolation of values DZ into the matrix P");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);

        //3 Tensioning and smoothing of the matrix P
        tension_loop(&mut abos_mutable, &abos_immutable);
        //println!("--------//3 tensioning");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);

        linear_tension_loop(&mut abos_mutable, &abos_immutable);
        //println!("--------//3 linear tensioning");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);

        smoothing_loop(&mut abos_mutable, &abos_immutable);
        //println!("--------//3 smoothing");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);

        //4 p + dp -> P
        abos_mutable.p += &abos_mutable.dp;
        //println!("--------//4 p + dp -> P");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);
        //back to 4

        //5 Dz - p => dz (Z - f (Xi, Yi) →DZi )
        calculate_dz(&mut abos_mutable, &abos_immutable);
        //println!("--------//5 Dz - p => dz (Z - f (Xi, Yi) →DZi )");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);
        // abos_mutable.dz -= abos_mutable.p;

        //6 calculate maximal difference
        let max_difference = abos_mutable.dz.max();
        println!("max_difference {}", max_difference);
        if max_difference.abs() < 0.1 {
            break;
        }
        //println!("--------//6 calculate maximal difference {}", max_difference);
        //output_all_matrixes(&&abos_mutable, &abos_immutable);

        //7 p -> Dp
        abos_mutable.dp.copy_from(&abos_mutable.p);
        //println!("--------//7 p -> Dp");
        //output_all_matrixes(&&abos_mutable, &abos_immutable);
    }
    output_all_matrixes(&&abos_mutable, &abos_immutable);
}

fn calculate_dz(abos_mutable: &mut ABOSMutable, abos_immutable: &ABOSImmutable) {
    // Uses bilinear interpolation
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

//
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

//makes a weighted average of 4 corner cells, top right, top left, bot right, bot left
//TODO shift from corner check, to mid perimeter check
fn tension_cell(ii: i32, jj: i32, k_i_j_mod: i32, abos_mutable: &mut ABOSMutable) -> f64 {
    //we need to get Pi , j=Pik , jPi , jkPi−k , jPi , j−k
    let p = &abos_mutable.p;
    let min_i: i32 = ii - k_i_j_mod;
    let min_j: i32 = jj - k_i_j_mod;
    let max_i: i32 = ii + k_i_j_mod;
    let max_j: i32 = jj + k_i_j_mod;

    let mut p1: f64 = 0.0;
    let mut p_counter = 0;
    unsafe {
        if min_i >= 0 && min_j >= 0 {
            p1 += *p.get_unchecked((min_i as usize, min_j as usize));
            p_counter += 1;
        }
        if max_i < p.nrows() as i32 && (max_j as usize) < p.ncols() {
            p1 += *p.get_unchecked((max_i as usize, max_j as usize));
            p_counter += 1;
        }
        if min_i >= 0 && max_j < p.ncols() as i32 {
            p1 += *p.get_unchecked((min_i as usize, max_j as usize));
            p_counter += 1;
        }
        if max_i < p.nrows() as i32 && min_j >= 0 {
            p1 += *p.get_unchecked((max_i as usize, min_j as usize));
            p_counter += 1;
        }
        if p_counter == 0 {
            p1 = -INFINITY;
        } else {
            p1 /= p_counter as f64;
        }
    }
    p1
}

pub fn per_parts_constant_interpolation(
    abos_mutable: &mut ABOSMutable,
    abos_immutable: &ABOSImmutable,
) {
    for (ii, row) in abos_immutable.nb.row_iter().enumerate() {
        for (jj, col) in row.iter().enumerate() {
            let point_closest = abos_immutable.xyz_points.row(*col);
            unsafe {
                let p_position = abos_mutable.p.get_unchecked_mut((ii, jj));
                *p_position = point_closest[2];
            }
        }
    }
}

pub fn tension_loop(abos_mutable: &mut ABOSMutable, abos_immutable: &ABOSImmutable) {
    let n = cmp::max(4, abos_immutable.k_max / 2 + 2);
    for n_countdown in (1..n + 1).rev() {
        for (ii, row) in abos_immutable.k_u_v.row_iter().enumerate() {
            for (jj, col) in row.iter().enumerate() {
                let k_i_j_mod = cmp::min(col.0, n_countdown);
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

fn linear_tension_cell(
    q: f64,
    ii: usize,
    jj: usize,
    u: usize,
    v: usize,
    abos_immutable: &ABOSImmutable,
    abos_mutable: &mut ABOSMutable,
) -> f64 {
    // println!("q {} u  {} v {}", q, u, v);
    let p = &abos_mutable.p;
    let top_l = if (ii + u) >= p.nrows() as usize || (jj + v) >= p.ncols() {
        0.0
    } else {
        //unsafe { *p.get_unchecked((ii + u, jj + v)) }
        unsafe { *abos_mutable.p.get_unchecked((ii + u, jj + v)) }
    };

    let bot_r = if (ii as i32 - u as i32) < 0 || (jj as i32 - v as i32) < 0 {
        0.0
    } else {
        //unsafe { *p.get_unchecked((ii - u, jj - v)) }
        unsafe { *abos_mutable.p.get_unchecked((ii - u, jj - v)) }
    };

    let top_r = if (ii as i32 - u as i32) < 0 || jj + v >= p.ncols() {
        0.0
    } else {
        //unsafe { *p.get_unchecked((ii - u, jj + v))}
        unsafe { *abos_mutable.p.get_unchecked((ii - u, jj + v)) }
    };

    let bot_l = if ii + u >= p.nrows() as usize || (jj as i32 - v as i32) < 0 {
        0.0
    } else {
        //unsafe { *p.get_unchecked((ii + u, jj - v))}
        unsafe { *abos_mutable.p.get_unchecked((ii + u, jj - v)) }
    };
    //code was
    //let bottom_divider = (q * top_l_b) + (q * top_r_b) + bot_l_b + bot_r_b;
    let bottom_divider = 2.0 * q + 2.0 * abos_immutable.r as f64;
    if bottom_divider == 0.0 {
        -INFINITY
    } else {
        (q * top_l + q * bot_r + top_r * abos_immutable.r as f64 + bot_l * abos_immutable.r as f64)
            / (bottom_divider)
    }
}

pub fn linear_tension_loop(abos_mutable: &mut ABOSMutable, abos_immutable: &ABOSImmutable) {
    let n = cmp::max(4, abos_immutable.k_max / 2 + 2);
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

fn get_q_linear_tension(abos_immutable: &ABOSImmutable, k_i_j: usize) -> f64 {
    //println!("self.l {}, self.k_max {}", self.l, self.k_max);
    match abos_immutable.degree {
        3 => 1.0,
        2 => abos_immutable.l * (abos_immutable.k_max - k_i_j) as f64,
        _ => abos_immutable.l * ((abos_immutable.k_max - k_i_j) as f64).powi(2),
    }
}

//tt : midpoint, checking window from
//dt : amount to look above below
//min_i : minimum acceptable value, inclusive
//max_i  : maximum acceptable value exclusive
// returns (tt - dt, tt + dt) pending on min_t, max_t
pub fn get_valid_dim_bounds(tt: usize, dt: usize, min_t: usize, max_t: usize) -> (usize, usize) {
    let lower_bound = if tt < dt || tt - dt < min_t {
        min_t
    } else {
        tt - dt
    };
    let upper_bound = if tt + dt > max_t { max_t } else { tt + dt };
    (lower_bound, upper_bound)
}

pub fn set_t_smooth(abos_mutable: &mut ABOSMutable) {
    //Set initial tt matrix
    for (ii, row) in abos_mutable.p.row_iter().enumerate() {
        let (kk_min, kk_max) = get_valid_dim_bounds(ii, 2, 0, abos_mutable.p.nrows());
        for (jj, _col) in row.iter().enumerate() {
            let (ll_min, ll_max) = get_valid_dim_bounds(jj, 2, 0, abos_mutable.p.ncols());
            //println!("kk_max {} kk_min {} ll_max {} ll_min {} ", kk_max, kk_min, ll_max, ll_min);
            let num_cells = ((kk_max - kk_min + 1) * (ll_max - ll_min + 1)) as f64;

            //.slice(start, shape)
            let mut pij_resid_sum = abos_mutable
                .p
                .slice((kk_min, ll_min), (kk_max - kk_min, ll_max - ll_min))
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
    abos_mutable.t_smooth.apply(|x| (x - min_t) * 1.0 / dt);
}

pub fn smoothing_loop(abos_mutable: &mut ABOSMutable, abos_immutable: &ABOSImmutable) {
    let n = cmp::max(4, abos_immutable.k_max * (abos_immutable.k_max / 16));

    for n_countdown in (1..n + 1).rev() {
        //t , are weights, which are zero before the first smoothing and afterwards they are computed according to the formula
        // abos_mutable.t_smooth.fill(0.0);
        set_t_smooth(abos_mutable);
        for (ii, row) in abos_immutable.k_u_v.row_iter().enumerate() {
            for (jj, col) in row.iter().enumerate() {
                let k_i_j = col.0;
                if k_i_j >= n_countdown {
                    let (kk_min, kk_max) =
                        get_valid_dim_bounds(ii, 1, 0, abos_mutable.p.nrows() - 1);
                    let (ll_min, ll_max) =
                        get_valid_dim_bounds(jj, 1, 0, abos_mutable.p.ncols() - 1);

                    let num_cells = ((kk_max - kk_min + 1) * (ll_max - ll_min + 1)) as f64;
                    let mut pij_new = abos_mutable
                        .p
                        .slice((kk_min, ll_min), (kk_max - kk_min, ll_max - ll_min))
                        .sum();
                    unsafe {
                        pij_new += num_cells
                            * abos_immutable.q_smooth
                            * (*abos_mutable.p.get_unchecked((ii, jj)))
                            * (*abos_mutable.t_smooth.get_unchecked((ii, jj)) - 1.0);
                        //num_cells - 1.0 inferred from case of equation 2.2.7
                        pij_new /= abos_immutable.q_smooth
                            * (*abos_mutable.t_smooth.get_unchecked((ii, jj)))
                            + (num_cells - 1.0);
                        *abos_mutable.p.get_unchecked_mut((ii, jj)) = pij_new;
                    }
                }
            }
        }
    }
}

pub fn output_all_matrixes(abos_mutable: &ABOSMutable, _abos_immutable: &ABOSImmutable) {
    //println!("dmc {}", abos_immutable.dmc);
    // println!(
    //     "NB {:.1} Z{:.1} DZ{:.1} DP{:.1} P{:.1}",
    //     abos_immutable.nb, abos_immutable.z, abos_mutable.dz, abos_mutable.dp, abos_mutable.p
    // );
    println!("{}", abos_mutable.p.len())
}

//
// #[cfg(test)]
// mod tests {
//     use crate::{ABOSGrid, swap_and_get_ranges, initialize_dmatrix,
//                 initialize_kdtree_from_matrix, get_min_chebyshev_distance_kdi};
//
//     extern crate nalgebra as na;
//
//     use na::{MatrixMN, Dynamic, U3};
//
//     extern crate rand;
//
//     use rand::prelude::*;
//
//
//     #[test]
//     fn test_swap_and_get_ranges() {
//         let points_fix = na::Matrix3::new(1.0, 2.0, 3.0,
//                                           2.0, 4.0, 6.0,
//                                           3.0, 6.0, 9.0);
//         let mut xyz_points: MatrixMN<f64, Dynamic, U3> = na::convert(points_fix);
//         let (x1, x2, y1, y2, z1, z2, xy_swaped) = swap_and_get_ranges(&mut xyz_points);
//
//         let check_points_fix = na::Matrix3::new(2.0, 1.0, 3.0,
//                                                 4.0, 2.0, 6.0,
//                                                 6.0, 3.0, 9.0);
//         let check_points: MatrixMN<f64, Dynamic, U3> = na::convert(check_points_fix);
//
//         assert_eq!(xyz_points, check_points);
//         assert_eq!((x1, x2, y1, y2, z1, z2, xy_swaped), (2.0, 6.0, 1.0, 3.0, 3.0, 9.0, true));
//     }
//
//     #[test]
//     fn test_initialize_dmatrix() {
//         //making initial 2d vector
//         let mut points: Vec<Vec<f64>> = vec![];
//         for ii in 0..3 {
//             //making f64 then converting seemse better than casting f64 3 times
//             let iif = 1.0 + ii as f64;
//             let new_point: Vec<f64> = vec!(iif, iif * 2.0, iif * 3.0);
//             //let new_point:Vec<f64> = vec!(ii , (ii*2) as f64, (ii*3) as f64);
//             points.push(new_point);
//         }
//
//         let xyz_points: MatrixMN<f64, Dynamic, U3> = initialize_dmatrix(points);
//
//         let check_xyz = na::Matrix3::new(1.0, 2.0, 3.0,
//                                          2.0, 4.0, 6.0,
//                                          3.0, 6.0, 9.0);
//         let check_xyz: MatrixMN<f64, Dynamic, U3> = na::convert(check_xyz);
//
//         //step 1: make an array with all the points
//         assert_eq!(xyz_points, check_xyz);
//     }
//
//     #[test]
//     fn test_min_chebyshev_dist() {
//         let check_xyz = na::Matrix3::new(1.0, 2.0, 3.0,
//                                          2.0, 4.0, 6.0,
//                                          3.0, 6.0, 9.0);
//         let check_xyz: MatrixMN<f64, Dynamic, U3> = na::convert(check_xyz);
//         let kdtree = initialize_kdtree_from_matrix(&check_xyz);
//         let cheby_dist = get_min_chebyshev_distance_kdi(&check_xyz, &kdtree);
//         assert_eq!(2.0, cheby_dist);
//     }
// }
