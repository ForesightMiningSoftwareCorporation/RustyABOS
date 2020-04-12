///Crate to generate a surface with the ABOS algorithm proposed in Art of Surface Interpolation byM gr. Miroslav Dressler, Ph.D.
/// http://m.dressler.sweb.cz/ABOS.htm
/// 
extern crate approx;
extern crate nalgebra as na;

pub const INFINITY: f64 = 1.0f64 / 0.0f64;

use kdtree::KdTree;
use kdtree::distance::squared_euclidean;
//use alga::general::{Real}; //generic index functions
use na::{DMatrix, DVector, MatrixMN, Dynamic, RealField, U3, U1, Dim};
mod deprecated;

//Make D matrix
//unwrapping to 1d vector [x1, y1, z1, x2, y2, z2...]
//real function call would just pass the vector into DMatrix
fn initialize_dmatrix(points: Vec<Vec<f64>>) -> MatrixMN<f64, Dynamic, U3> {
    let mut vec = Vec::new();
    for point in points.iter() {
        for ii in point.iter() {
            vec.push(*ii);
        }
    }
    let point_count = points.len();
    let dm: DMatrix<f64> = DMatrix::from_iterator(3, point_count, vec.into_iter());
    let fix_dm = dm.transpose();

    //step 1: make an array with all the points
    let xyz_points: MatrixMN<f64, Dynamic, U3> = na::convert(fix_dm);
    xyz_points
}

fn initialize_kdtree_from_matrix(xyz_points: &MatrixMN<f64, Dynamic, U3>) -> KdTree<f64, usize, [f64; 2]> {
    let mut kdtree = KdTree::new(2);
    // let mut points:Vec<([f64; 2], usize)> = vec![];
    for (i, row) in xyz_points.row_iter().enumerate() {
        kdtree.add([row[0], row[1]], i).unwrap();
    }
    kdtree
}

fn get_min_chebyshev_distance_kdi(xyz_points: &MatrixMN<f64, Dynamic, U3>, kdtree : &KdTree<f64, usize, [f64; 2]>) -> f64 {
    let mut min_chebyshev_distance: f64 = INFINITY;

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

fn get_scaled_u_v(u : f64, v : f64, n : f64) -> (f64, f64){
    let (mut u_mod, mut v_mod) = (u, v);
    let uv_magnitude = (u*u + v*v).sqrt();
    if uv_magnitude > n {
        let c = n / uv_magnitude;
        u_mod = c * u;
        v_mod = c * v;
    };
    (u_mod, v_mod)
}

//makes a weighted average of 4 corner cells, top right, top left, bot right, bot left
fn linear_tension_cell(q : f64, ii : usize, jj : usize, u : usize, v :usize, p: &mut MatrixMN<f64, Dynamic, Dynamic>) -> f64 {
    //println!("q {} u  {} v {}", q, u, v);
    let top_l = if ii + u >= p.nrows() as usize || jj + v >= p.ncols() { 
        0.0
    } else {unsafe { *p.get_unchecked_mut((ii + u, jj + v)) } };

    let bot_r = if (ii as i32 - u as i32) < 0  || (jj as i32 - v as i32) < 0 { 
        0.0
    } else { unsafe {*p.get_unchecked_mut((ii - u, jj - v)) }};

    let top_r = if (ii as i32 - u as i32) < 0 || jj + v >= p.ncols() { 
        0.0
    } else {unsafe { *p.get_unchecked_mut((ii - u, jj + v)) } };

    let bot_l = if ii + u >= p.nrows() as usize  || (jj as i32 - v as i32) < 0 { 
        0.0
    } else { unsafe {*p.get_unchecked_mut((ii + u, jj - v)) }};
    //println!("top_l {} bot_r {} top_r {} bot_l {}", top_l, bot_r, top_r, bot_l);
    let lin_ten = (q*(top_l * bot_r) + top_r + bot_l)/(2.0 * q + 2.0);
    //println!("lin_ten {} ", lin_ten);
    lin_ten
}



fn tension_cell(ii: i32, jj: i32, k_i_j_mod: i32, p: &mut MatrixMN<f64, Dynamic, Dynamic>) -> f64 {
    //we need to get Pi , j=Pik , jPi , jkPi−k , jPi , j−k
    let min_i: i32 = ii - k_i_j_mod;
    let min_j: i32 = jj - k_i_j_mod;
    let max_i: i32 = ii + k_i_j_mod;
    let max_j: i32 = jj + k_i_j_mod;


    let mut p1: f64 = 0.0;
    let mut p_counter = 0;
    unsafe {
        if min_i >= 0 && min_j >= 0 {
            p1 += *p.get_unchecked((min_i as usize, min_j as usize));
            p_counter+=1;
        }
        if (max_i as usize) < p.nrows() && (max_j as usize) < p.ncols() {
            p1 += *p.get_unchecked((max_i as usize, max_j as usize));
            p_counter+=1;
        }
        if min_i >= 0 && (max_j as usize) < p.ncols() {
            p1 += *p.get_unchecked((min_i as usize, max_j as usize));
            p_counter+=1;
        }
        if (max_i as usize) < p.nrows() && min_j >= 0 {
            p1 += *p.get_unchecked((max_i as usize, min_j as usize));
            p_counter+=1;
        }

        p1 = p1 / p_counter as f64;
    }
    p1
}


pub struct ABOSGrid {
    //User Inputs
    degree: i8,
    r: usize,
    l: f64,
    filter: f64,
    pub n: usize,
    //INPUT resolution parameter
    xyz_points: MatrixMN<f64, Dynamic, U3>,
    //INPUT all points XYZ

    //Derived Points
    x1: f64,
    //minx
    x2: f64,
    //maxx
    y1: f64,
    //miny
    y2: f64,
    //maxy
    z1: f64,
    //min z
    z2: f64,
    //max z

    //ABOS Calculated Parameters
    dmc: f64,
    //minimal chebyshev distance
    pub i1: i32,
    //xsize of grid
    pub j1: i32,
    //ysize of grid
    dx: f64,
    //size of grid on x
    dy: f64,
    //size of grid on y
    pub p: MatrixMN<f64, Dynamic, Dynamic>,
    //elements of the elevation matrix [i,j] size
    dp: MatrixMN<f64, Dynamic, Dynamic>,
    //auxiliary matrix same size as P
    nb: MatrixMN<usize, Dynamic, Dynamic>,
    //[i x j]
    z: DVector<f64>,
    //vector if z coordinates XYZ
    dz: DVector<f64>,
    //auxiliary vector same size as Z
    pub k_u_v: MatrixMN<(usize, usize, usize), Dynamic, Dynamic>, //first one is max, second one is x dist, third one is y dist
    // Grid distance of each grid to the point indexed in NB
    k_max: usize,
    //maximal element of matrix K
    rs: f64,
    // Resolution of map
    xy_swaped: bool,
    //whether the xy points were swapped
}

impl ABOSGrid {
    //points 2d vectory of [[x1,y1,z2], [x2,y2,z2]....]
    pub fn new(points: Vec<Vec<f64>>, filter: f64, degree: i8) -> ABOSGrid {
        //unwrapping to 1d vector [x1, y1, z1, x2, y2, z2...]
        //real function call would just pass the vector into DMatrix
        let mut vec = Vec::new();
        for point in points.iter() {
            for ii in point.iter() {
                vec.push(*ii);
            }
        }
        let point_count = points.len();
        let dm: DMatrix<f64> = DMatrix::from_iterator(3, point_count, vec.into_iter());
        let fix_dm = dm.transpose();

        //step 1: make an array with all the points
        let mut xyz_points: MatrixMN<f64, Dynamic, U3> = initialize_dmatrix(points);

        //step 2: get Point range information and swap as necessary
        let (x1, x2, y1, y2, z1, z2, xy_swaped) = swap_and_get_ranges(&mut xyz_points);

        let kdtree = initialize_kdtree_from_matrix(&xyz_points);
        //step 3: get Chebyshev distance
        let dmc = get_min_chebyshev_distance_kdi(&xyz_points, &kdtree);
        //step 3: get the grid dimensions
        let (i1, j1, dx, dy) = compute_grid_dimensions(x1, x2, y1, y2, dmc, filter);

        //step 4: Create empty vectors
        let p: MatrixMN<f64, Dynamic, Dynamic> = MatrixMN::from_element_generic(Dynamic::from_usize(i1 as usize), Dynamic::from_usize(j1 as usize), 0.0);
        let dp: MatrixMN<f64, Dynamic, Dynamic> = p.clone_owned();
        let nb: MatrixMN<usize, Dynamic, Dynamic> = MatrixMN::from_element_generic(Dynamic::from_usize(i1 as usize), Dynamic::from_usize(j1 as usize), 0);
        let k_u_v: MatrixMN<(usize, usize, usize), Dynamic, Dynamic> = MatrixMN::from_element_generic(Dynamic::from_usize(i1 as usize), Dynamic::from_usize(j1 as usize), (0, 0 ,0));

        // Pcontainer.matrix::
        let z: DVector<f64> = xyz_points.column(2).clone_owned();
        let dz = z.clone_owned();

        let res_x = (x2 - x1) / filter;
        let res_y = (y2 - y1) / filter;
        let rs = if res_x > res_y { res_x } else { res_y };

        // These items must be calculated with calculated k_max
        let k_max = 0;
        let r = 0;
        let l = 0.0;
        let n = 0;
        let mut abos_grid = ABOSGrid {
            degree,
            r,
            l,
            n,
            filter, //INPUT resolution parameter
            xyz_points, //INPUT all points XYZ
            x1, //minx
            x2, //maxx
            y1, //miny
            y2, //maxy
            z1, //min z
            z2, //max z
            dmc, //minimal chebyshev distance
            i1, //xsize of grid
            j1, //ysize of grid
            dx, //size of grid on x
            dy, //size of grid on y
            p, //elements of the elevation matrix [i,j] size
            dp, //auxiliary matrix same size as p
            z, //vector if z coordinates XYZ
            dz, //auxiliary vector same size as z
            nb, // Matrix of nearest points on grid. Containing indexes to nearest point in the XYZ array
            k_u_v, // Grid distance of each grid to the point indexed in nb
            k_max,  //maximal element of matrix k
            rs, // Resolution of map
            xy_swaped,
        };
        //Calculates k_u_v and kmax
        abos_grid.init_distance_point_matrixes_kdi(&kdtree);

        let (r, l) = compute_rl(degree, abos_grid.k_max);
        abos_grid.r = r;
        abos_grid.l = l;
        abos_grid.n = std::cmp::max(4, k_max / 2 + 2);
        println!("k_max {} r {}, l {}, n {}",  abos_grid.k_max, r, l , n);
        abos_grid
    }

    pub fn per_parts_constant_interpolation(&mut self) {
        for (ii, row) in self.nb.row_iter().enumerate() {
            for (jj, col) in row.iter().enumerate() {
                let point_closest = self.xyz_points.row(*col);
                unsafe {
                    let p_position = self.p.get_unchecked_mut((ii, jj));
                    *p_position = point_closest[2];
                }
            }
        }
    }

    fn indexes_to_position(&self, row_index: &usize, col_index: &usize) -> [f64; 2] {
        let x_coordinate = self.x1 + self.dx * (*row_index as f64);
        let y_coordinate = self.y1 + self.dy * (*col_index as f64);

        [x_coordinate, y_coordinate]
    }


    pub fn tension_loop(abos_grid: &mut ABOSGrid) {
        for n_countdown in (1..abos_grid.n + 1).rev() {
            for (ii, row) in abos_grid.k_u_v.row_iter().enumerate() {
                for (jj, col) in row.iter().enumerate() {
                    let k_i_j_mod = std::cmp::min(col.0, n_countdown);
                    unsafe {
                        *abos_grid.p.get_unchecked_mut((ii, jj)) = tension_cell(ii as i32, jj as i32, k_i_j_mod as i32, &mut abos_grid.p);
                    }
                }
            }
        }
    }


    pub fn linear_tension_loop(abos_grid: &mut ABOSGrid) {
        for n_countdown in (1..abos_grid.n + 1).rev() {
            for (ii, row) in abos_grid.k_u_v.row_iter().enumerate() {
                for (jj, col) in row.iter().enumerate() {
                    let k_i_j_mod = std::cmp::min(col.0, n_countdown);
                    let q = abos_grid.get_q_value(k_i_j_mod);
                    //println!("q {}", q);
                    let(u_mod, v_mod) = get_scaled_u_v(col.1 as f64, col.2 as f64, abos_grid.n as f64);

                    unsafe {
                        *abos_grid.p.get_unchecked_mut((ii, jj)) = linear_tension_cell(q, u_mod as usize, v_mod as usize, ii, jj, &mut abos_grid.p);
                    }
                }
            }
        }
    }
    
    fn get_q_value(&self, k_i_j: usize) -> f64 {
        //println!("self.l {}, self.k_max {}", self.l, self.k_max);
        return match self.degree {
            3 => 1.0,
            2 => {
                self.l * (self.k_max - k_i_j) as f64
            }
            _ => {
                self.l * ((self.k_max - k_i_j) as f64).powi(2)
            }
        };
    }


    pub fn output_all_matrixes(&self) {
        println!("dmc {}", self.dmc);
        println!("NB {:.1} Z{:.1} DZ{:.1} DP{:.1} P{:.1}", self.nb, self.z, self.dz, self.dp, self.p);
    }


    pub fn init_distance_point_matrixes_kdi(&mut self, kdtree : &KdTree<f64, usize, [f64; 2]>) {
        //iterate through each grid cell position. Set index of point to NB, and grid distance to K
        for (ii, row) in self.dp.row_iter().enumerate() {
            for (jj, _col) in row.iter().enumerate() {
                let position = self.indexes_to_position(&ii, &jj);
                let kd_search_result = kdtree.nearest(&position, 1, &squared_euclidean).unwrap();
                let closest_point_in_tree_indx = *kd_search_result[0].1;

                let closest_point = self.xyz_points.row(closest_point_in_tree_indx);
                let x_distance = f64::round(f64::abs((closest_point[0] - position[0]) / self.dx)) as usize;
                let y_distance = f64::round(f64::abs((closest_point[1] - position[1]) / self.dy)) as usize;

                unsafe {
                    let nb_position = self.nb.get_unchecked_mut((ii, jj));
                    *nb_position = closest_point_in_tree_indx;

                    let k_position = self.k_u_v.get_unchecked_mut((ii, jj));
                    //*k_position = if x_distance > y_distance { x_distance } else { y_distance };
                    let max_cell_dist = if x_distance > y_distance { x_distance } else { y_distance };
                    *k_position = (max_cell_dist, x_distance, y_distance);
                    if k_position.0 > (self.k_max) {
                        self.k_max = k_position.0;
                    }
                }
            }
        }
    }

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
    
    pub fn calculation_loop(&mut self){
        let mut numLoops = self.n;
        //numLoops = 1;
        while numLoops > 0 {
            self.per_parts_constant_interpolation();
            ABOSGrid::tension_loop(self);
            ABOSGrid::linear_tension_loop(self);
            numLoops -= 1;
        }
    }

} //end abos


pub fn compute_grid_dimensions(x1: f64, x2: f64, y1: f64, y2: f64, dmc: f64, filter: f64)
                               -> (i32, i32, f64, f64) {
    //step 1: Always assuming x side is greater

    //step 2: grid size is defined as i0=round  x21/ Dmc 

    let i0: i32 = f64::round((x2 - x1) / dmc) as i32;
    //step 3 find your grid size for your larger dimension set as i1 = i0*k largest possible while less than Filter
    let mut i1: i32 = 0;
    let mut i = 0;
    loop {
        i += 1;
        let potential_val: i32 = i0 * i;
        if filter > potential_val as f64 {
            i1 = potential_val;
        } else {
            break;
        }
    }

    //step 5 find your grid size for your smaller dimension such that it is as close to that of il as possible
    let j1: i32 = f64::round((y2 - y1) / (x2 - x1) * (i1 as f64 - 1.0)) as i32;

    let dx = (x2 - x1) / (i1 - 1) as f64; //include the minus one so matices inclde the max points
    let dy = (y2 - y1) / (j1 - 1) as f64;

    return (i1, j1, dx, dy);
}

fn compute_rl(degree: i8, k_max: usize) -> (usize, f64) {
    //println!();
    return match degree {
        0 => {
            let r = 1;
            let l = 0.7 / ((0.107 * k_max as f64 - 0.714) * k_max as f64);
            (r, l)
        }
        1 => {
            let r = 1;
            let l = 1.0 / ((0.107 * k_max as f64 - 0.714) * k_max as f64);
            (r, l)
        }
        2 => {
            let r = 1;
            let l = 1.0 / (0.0360625 * k_max as f64 + 0.192);
            (r, l)
        }
        3 => {
            let r = 0;
            let l = 0.7 / ((0.107 * k_max as f64 - 0.714) * k_max as f64);
            (r, l)
        }
        _ => {
            (0, 0.0)
        }
    };
}

pub fn get_ranges(points: &MatrixMN<f64, Dynamic, U3>) -> (f64, f64, f64, f64, f64, f64) {
    let x1 = points.column(0).min();
    let x2 = points.column(0).max();
    let y1 = points.column(1).min();
    let y2 = points.column(1).max();
    let z1 = points.column(2).min();
    let z2 = points.column(2).max();

    return (x1, x2, y1, y2, z1, z2);
}

pub fn swap_and_get_ranges(points: &mut MatrixMN<f64, Dynamic, U3>) -> (f64, f64, f64, f64, f64, f64, bool) {
    let (x1, x2, y1, y2, z1, z2) = get_ranges(&points);
    return if f64::abs(x1 - x2) < f64::abs(y1 - y2) {
        points.swap_columns(0, 1);
        let (x1, x2, y1, y2, z1, z2) = get_ranges(&points);
        (x1, x2, y1, y2, z1, z2, true)
    } else {
        (x1, x2, y1, y2, z1, z2, false)
    };
}

#[cfg(test)]
mod tests {
    use crate::{ABOSGrid, swap_and_get_ranges, initialize_dmatrix, 
                initialize_kdtree_from_matrix, get_min_chebyshev_distance_kdi};

    extern crate nalgebra as na;

    use na::{MatrixMN, Dynamic, U3};
    extern crate rand;
    use rand::prelude::*;



    #[test]
    fn test_swap_and_get_ranges() {
        let points_fix = na::Matrix3::new(1.0, 2.0, 3.0,
                                          2.0, 4.0, 6.0,
                                          3.0, 6.0, 9.0);
        let mut xyz_points: MatrixMN<f64, Dynamic, U3> = na::convert(points_fix);
        let (x1, x2, y1, y2, z1, z2, xy_swaped) = swap_and_get_ranges(&mut xyz_points);

        let check_points_fix = na::Matrix3::new(2.0, 1.0, 3.0,
                                                4.0, 2.0, 6.0,
                                                6.0, 3.0, 9.0);
        let check_points: MatrixMN<f64, Dynamic, U3> = na::convert(check_points_fix);

        assert_eq!(xyz_points, check_points);
        assert_eq!((x1, x2, y1, y2, z1, z2, xy_swaped), (2.0, 6.0, 1.0, 3.0, 3.0, 9.0, true));
    }

    #[test]
    fn test_initialize_dmatrix() {
        //making initial 2d vector
        let mut points: Vec<Vec<f64>> = vec![];
        for ii in 0..3 {
            //making f64 then converting seemse better than casting f64 3 times
            let iif = 1.0 + ii as f64;
            let new_point: Vec<f64> = vec!(iif, iif * 2.0, iif * 3.0);
            //let new_point:Vec<f64> = vec!(ii , (ii*2) as f64, (ii*3) as f64);
            points.push(new_point);
        }

        let xyz_points: MatrixMN<f64, Dynamic, U3> = initialize_dmatrix(points);

        let check_xyz = na::Matrix3::new(1.0, 2.0, 3.0,
                                         2.0, 4.0, 6.0,
                                         3.0, 6.0, 9.0);
        let check_xyz: MatrixMN<f64, Dynamic, U3> = na::convert(check_xyz);

        //step 1: make an array with all the points
        assert_eq!(xyz_points, check_xyz);
    }

    #[test]
    fn test_min_chebyshev_dist() {

        let check_xyz = na::Matrix3::new(1.0, 2.0, 3.0,
                                         2.0, 4.0, 6.0,
                                         3.0, 6.0, 9.0);
        let check_xyz: MatrixMN<f64, Dynamic, U3> = na::convert(check_xyz);
        let kdtree = initialize_kdtree_from_matrix(&check_xyz);
        let cheby_dist = get_min_chebyshev_distance_kdi(&check_xyz, &kdtree);
        assert_eq!(2.0, cheby_dist);
    }
}