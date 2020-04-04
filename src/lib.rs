///Crate to generate a surface with the ABOS algorithm proposed in Art of Surface Interpolation byM gr. Miroslav Dressler, Ph.D.
/// http://m.dressler.sweb.cz/ABOS.htm
/// 

extern crate approx;
extern crate nalgebra as na;

pub const INFINITY: f64 = 1.0f64 / 0.0f64;

use kdtree::KdTree;
use kdtree::ErrorKind;
use kdtree::distance::squared_euclidean;

use na::{DMatrix, DVector, MatrixMN, Dynamic, U3, U1, Dim};

//TODO make closure 
fn get_min_chebyshev_distance(xyz_points: &MatrixMN<f64, Dynamic, U3>) -> f64 {
    let mut min_chebyshev_distance: f64 = INFINITY;
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
    degree: i8,
    r: usize,
    l: f64,
    filter: f64,
    n: usize,
    //INPUT resolution parameter
    xyz_points: MatrixMN<f64, Dynamic, U3>,
    //INPUT all points XYZ
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
    dmc: f64,
    //minimal chebyshev distance
    i1: i32,
    //xsize of grid
    j1: i32,
    //ysize of grid
    dx: f64,
    //size of grid on x
    dy: f64,
    //size of grid on y
    p: MatrixMN<f64, Dynamic, Dynamic>,
    //elements of the elevation matrix [i,j] size
    dp: MatrixMN<f64, Dynamic, Dynamic>,
    //auxiliary matrix same size as P
    nb: MatrixMN<usize, Dynamic, Dynamic>,
    //[i x j]
    z: DVector<f64>,
    //vector if z coordinates XYZ
    dz: DVector<f64>,
    //auxiliary vector same size as Z
    k: MatrixMN<usize, Dynamic, Dynamic>,
    // Grid distance of each grid to the point indexed in NB
    k_max: usize,
    //maximal element of matrix K
    rs: f64, // Resolution of map
}

impl ABOSGrid {
    pub fn new(points: Vec<Vec<f64>>, filter: f64, degree: i8) -> ABOSGrid {
        let mut vec = Vec::new(); //real function call would just pass the vector in here
        for x in points.iter() {
            for y in x.iter() {
                vec.push(*y);
            }
        }
        let point_count = points.len();
        let dm: DMatrix<f64> = DMatrix::from_iterator(3, point_count, vec.into_iter());
        let fix_dm = dm.transpose();

        //step 1: make an array with all the points
        let xyz_points: MatrixMN<f64, Dynamic, U3> = na::convert(fix_dm);
        //step 2: get Point range information
        let (x1, x2, y1, y2, z1, z2) = get_ranges(&xyz_points);

        //step 3: get Chebyshev distance
        let dmc = get_min_chebyshev_distance(&xyz_points);
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
        }
    }

    fn indexes_to_position(&self, row_index: &usize, col_index: &usize) -> [f64; 2] {
        let x_coordinate = self.x1 + self.dx * (*row_index as f64);
        let y_coordinate = self.y1 + self.dy * (*col_index as f64);

        [x_coordinate, y_coordinate]
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

    pub fn tension_loop(&mut self) {
        for n_countdown in (1..self.n + 1).rev() {
            self.tension_grid(&n_countdown);
        }
    }

    fn tension_grid(&mut self, n_countdown: &usize) {
        for (rowIndex, row) in self.k.row_iter().enumerate() {
            for (colIndex, col) in row.iter().enumerate() {
                let k_to_use = std::cmp::min(col,n_countdown);
            }
        }
    }

    fn tension_cell(rowIndex: &usize, colIndex: &usize, k_i_j_mod: &usize, p: &MatrixMN<f64, Dynamic, Dynamic>) {}

    pub fn output_all_matrixes(&self) {
        println!("NB {} K{} Z{} DZ{} DP{} P{:.1}", self.nb, self.k, self.z, self.dz, self.dp, self.p);
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
                let closest_point_in_tree = *kd_search_result[0].1;

                let closest_point = self.xyz_points.row(closest_point_in_tree);
                let x_distance = f64::round(f64::abs((closest_point[0] - position[0]) / self.dx)) as usize;
                let y_distance = f64::round(f64::abs((closest_point[1] - position[1]) / self.dy)) as usize;

                unsafe {
                    let nb_position = self.nb.get_unchecked_mut((rowIndex, colIndex));
                    *nb_position = closest_point_in_tree;

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

#[cfg(test)]
mod tests {
    use crate::{compute_grid_dimensions, ABOSGrid};

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
}


pub fn hello1() -> () {
    println!("Hello1");
}

pub fn hello2() -> () {
    println!("Hello2");
}

