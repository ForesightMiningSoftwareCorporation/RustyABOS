///Crate to generate a surface with the ABOS algorithm proposed in Art of Surface Interpolation byM gr. Miroslav Dressler, Ph.D.
/// http://m.dressler.sweb.cz/ABOS.htm
/// 

extern crate approx;
extern crate nalgebra as na;
//extern crate alga;

pub const INFINITY: f64 = 1.0f64 / 0.0f64;

use kdtree::KdTree;
use kdtree::ErrorKind;
use kdtree::distance::squared_euclidean;
//use alga::general::{Real}; //generic index functions
use na::{DMatrix, DVector, MatrixMN, Dynamic, RealField, U3, U1, Dim};

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

fn indexes_to_position_direct<N: RealField>(indxs: DVector<N>, d_mins: DVector<N>, d_sizes: DVector<N>) -> DVector<N> {
    let position = d_mins + d_sizes * indxs;
    position
}

fn indexes_to_position(indxs: DVector<u32>, d_mins: DVector<f64>, d_sizes: DVector<f64>) -> DVector<f64> {
    let indxsConverted: DVector<f64> = na::convert(indxs);
    indexes_to_position_direct(indxsConverted, d_mins, d_sizes)
}


//TODO convert from n2 to kd tree implimentation
fn get_min_chebyshev_distance_n2(xyz_points: &MatrixMN<f64, Dynamic, U3>) -> f64 {
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

fn get_min_cheyshev_distance_kd(xyz_points: &MatrixMN<f64, Dynamic, U3>) -> f64 {
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
        let kd_search_result = kdtree.nearest(&point, 1, &squared_euclidean).unwrap();
        let closest_indx = *kd_search_result[0].1;
        let distances: MatrixMN<f64, U1, U3> = row.clone_owned() - xyz_points.row(closest_indx);
        let distances = distances.abs();
        let max_xy_distance = if distances[0] > distances[1] { distances[0] } else { distances[1] };

        if max_xy_distance < min_chebyshev_distance {
            min_chebyshev_distance = max_xy_distance;
        }
        // unsafe {
        //     let distances = [point[0] - xyz_points[(closest_indx, 0)] ]
        //     *pPosition = point_closest[2];

        // }


        // for (colIndex, mut col) in row.iter().enumerate() {
        //     // let position = self.indexes_to_position(&rowIndex, &colIndex);
        //     // let kd_search_result = kdtree.nearest(&position, 1, &squared_euclidean).unwrap();
        //     // let closest_point_in_tree_indx = *kd_search_result[0].1;
        //     // let closest_point = self.xyz_points.row(closest_point_in_tree_indx);
        //     // let x_distance = f64::round(f64::abs((closest_point[0] - position[0]) / self.dx));
        //     // let y_distance = f64::round(f64::abs((closest_point[1] - position[1]) / self.dy));
        // }
    }

    let xyz_iter1 = xyz_points.row_iter();
    for (_, row) in xyz_iter1.enumerate() {
        let xyz_iter2 = xyz_points.row_iter();
        for (_, row2) in xyz_iter2.enumerate() {
            if row != row2 {
                let distances: MatrixMN<f64, U1, U3> = row.clone_owned() - row2;
                let distances = distances.abs();
                let max_local_distance = if distances[0] > distances[1] { distances[0] } else { distances[1] };

                if max_local_distance < min_chebyshev_distance {
                    min_chebyshev_distance = max_local_distance;
                }
            }
        }
    }
    return min_chebyshev_distance;
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
    pub k: MatrixMN<usize, Dynamic, Dynamic>,
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
        let mut xyz_points: MatrixMN<f64, Dynamic, U3> = na::convert(fix_dm);

        let mut xyz_points: MatrixMN<f64, Dynamic, U3> = initialize_dmatrix(points);

        //step 2: get Point range information and swap as necessary
        let (x1, x2, y1, y2, z1, z2, xy_swaped) = swap_and_get_ranges(&mut xyz_points);

        //step 3: get Chebyshev distance
        let dmc = get_min_chebyshev_distance_n2(&xyz_points);
        //step 3: get the grid dimensions
        let (i1, j1, dx, dy) = compute_grid_dimensions(x1, x2, y1, y2, dmc, filter);

        //step 4: Create empty vectors
        let p: MatrixMN<f64, Dynamic, Dynamic> = MatrixMN::from_element_generic(Dynamic::from_usize(i1 as usize), Dynamic::from_usize(j1 as usize), 0.0);
        let dp: MatrixMN<f64, Dynamic, Dynamic> = p.clone_owned();
        let nb: MatrixMN<usize, Dynamic, Dynamic> = MatrixMN::from_element_generic(Dynamic::from_usize(i1 as usize), Dynamic::from_usize(j1 as usize), 0);
        let k: MatrixMN<usize, Dynamic, Dynamic> = nb.clone_owned();

        // Pcontainer.matrix::
        let z: DVector<f64> = xyz_points.column(2).clone_owned();
        let dz = z.clone_owned();

        let k_max = 0;
        let res_x = (x2 - x1) / filter;
        let res_y = (y2 - y1) / filter;
        let rs = if res_x > res_y { res_x } else { res_y };
        // step 5: compute R and L
        let (r, l) = compute_rl(degree, k_max);
        let n = std::cmp::max(4, k_max / 2 + 2);
        ABOSGrid {
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
            k, // Grid distance of each grid to the point indexed in nb
            k_max,  //maximal element of matrix k
            rs, // Resolution of map
            xy_swaped,
        }
    }

    pub fn per_parts_constant_interpolation(&mut self) {
        for (rowIndex, row) in self.nb.row_iter().enumerate() {
            for (colIndex, col) in row.iter().enumerate() {
                let point_closest = self.xyz_points.row(*col);
                unsafe {
                    let pPosition = self.p.get_unchecked_mut((rowIndex, colIndex));
                    *pPosition = point_closest[2];
                }
            }
        }
    }

    fn get_q_value(&self, k_i_j: usize) -> f64 {
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

    fn indexes_to_position(&self, row_index: &usize, col_index: &usize) -> [f64; 2] {
        let x_coordinate = self.x1 + self.dx * (*row_index as f64);
        let y_coordinate = self.y1 + self.dy * (*col_index as f64);

        [x_coordinate, y_coordinate]
    }

    pub fn tension_loop( n: usize,i1:usize,j1:usize, mutable_p: &mut MatrixMN<f64, Dynamic, Dynamic>, k: &MatrixMN<usize, Dynamic, Dynamic>) {
        for n_countdown in (1..n + 1).rev() {
            for (rowIndex, row) in k.row_iter().enumerate() {
                for (colIndex, col) in row.iter().enumerate() {
                    let k_to_use = std::cmp::min(col, &n_countdown);
                    let new_p = ABOSGrid::tension_cell(i1,j1, rowIndex, colIndex, *k_to_use, mutable_p);
                    unsafe {
                        *mutable_p.get_unchecked_mut((rowIndex, colIndex)) = new_p;
                    }
                }
            }
        }
    }

    fn tension_cell(i1: usize,j1: usize, rowIndex: usize, colIndex: usize, k_i_j_mod: usize, p: &mut MatrixMN<f64, Dynamic, Dynamic>) -> f64 {
        //we need to get Pi , j=Pik , jPi , jkPi−k , jPi , j−k
        let min_i: usize = std::cmp::max((rowIndex as i32 - k_i_j_mod as i32), 0) as usize;
        let min_j: usize = std::cmp::max((colIndex as i32 - k_i_j_mod as i32), 0) as usize;
        let max_i: usize = std::cmp::min((rowIndex as i32 + k_i_j_mod as i32), (i1 as i32 - 1)) as usize;
        let max_j: usize = std::cmp::min((colIndex as i32 + k_i_j_mod as i32), (j1 as i32 - 1)) as usize;



        let mut p1: f64 = 0.0;
        unsafe {
            p1 += *p.get_unchecked((min_i, min_j));
            p1 += *p.get_unchecked((max_i, max_j));
            p1 += *p.get_unchecked((min_i, max_j));
            p1 += *p.get_unchecked((max_i, min_j));
            //
            p1 = p1 / 4.0;
        }
        p1
    }

    pub fn output_all_matrixes(&self) {
        println!("dmc {}", self.dmc);
        println!("NB {:.1} K{:.1} Z{:.1} DZ{:.1} DP{:.1} P{:.1}", self.nb, self.k, self.z, self.dz, self.dp, self.p);
    }

    pub fn init_distance_point_matrixes(&mut self) {
        //step 1: make a 2d search tree to accelerate the finding of points
        let dimensions = 2;
        let mut kdtree = KdTree::new(dimensions);
        // let mut points:Vec<([f64; 2], usize)> = vec![];
        for (i, row) in self.xyz_points.row_iter().enumerate() {
            kdtree.add([row[0], row[1]], i).unwrap();
        }
        //step 2: iterate through each grid cell position. Set index of point to NB, and grid distance to K
        for (rowIndex, mut row) in self.dp.row_iter().enumerate() {
            for (colIndex, mut col) in row.iter().enumerate() {
                let position = self.indexes_to_position(&rowIndex, &colIndex);
                let kd_search_result = kdtree.nearest(&position, 1, &squared_euclidean).unwrap();
                let closest_point_in_tree_indx = *kd_search_result[0].1;

                let closest_point = self.xyz_points.row(closest_point_in_tree_indx);
                let x_distance = f64::round(f64::abs((closest_point[0] - position[0]) / self.dx)) as usize;
                let y_distance = f64::round(f64::abs((closest_point[1] - position[1]) / self.dy)) as usize;

                unsafe {
                    let nb_position = self.nb.get_unchecked_mut((rowIndex, colIndex));
                    *nb_position = closest_point_in_tree_indx;

                    let k_position = self.k.get_unchecked_mut((rowIndex, colIndex));
                    *k_position = if x_distance > y_distance { x_distance } else { y_distance };
                    if *k_position > (self.k_max) {
                        self.k_max = *k_position
                    }
                }
            }
        }
    }
}


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

    let dx = (x2 - x1) / i1 as f64;
    let dy = (y2 - y1) / j1 as f64;

    return (i1, j1, dx, dy);
}

fn compute_rl(degree: i8, k_max: usize) -> (usize, f64) {
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
    use crate::{compute_grid_dimensions, ABOSGrid, swap_and_get_ranges, initialize_dmatrix, get_min_cheyshev_distance_kd};

    extern crate nalgebra as na;

    use na::{DMatrix, DVector, MatrixMN, Dynamic, U3, U1, Dim};
// #[test]
// fn it_works() {
//     assert_eq!(2 + 2, 4);
// }

// #[test]
// fn nalgebraImportCheck() {
//     let axis  = na::Vector3::x_axis();
//     let angle = 1.57;
//     let b     = na::Rotation3::from_axis_angle(&axis, angle);
//
//     approx::relative_eq!(b.axis().unwrap(), axis);
//     approx::relative_eq!(b.angle(), angle);
// }

    #[test]
    fn test_correct_grid_size() {
// let mut points: Vec<na::Vector3<f64>> = vec![];
//
// for i in 0..10 {
//     let new_point = na::Vector3::new(rng.gen(), rng.gen(), rng.gen());
//     points.push(new_point);
// }
// println!("{:?}", points);
// test_abos = ABOSGrid::new(points);
    }

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
        //let min_cheby_dist = get_min_chebyshev_distance_n2(&check_xyz);
        let cheby_dist = get_min_cheyshev_distance_kd(&check_xyz);
        println!("min_cheby_dist {}", get_min_cheyshev_distance_kd(&check_xyz));
        assert_eq!(2.0, cheby_dist);
    }
}


pub fn hello1() -> () {
    println!("Hello1");
}

pub fn hello2() -> () {
    println!("Hello2");
}

