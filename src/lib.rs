///Crate to generate a surface with the ABOS algorithm proposed in Art of Surface Interpolation byM gr. Miroslav Dressler, Ph.D.
/// http://m.dressler.sweb.cz/ABOS.htm
/// 

extern crate approx;
// For the macro relative_eq! in nalgebraImportCheck
extern crate nalgebra as na;

pub const INFINITY: f64 = 1.0f64 / 0.0f64;

use na::{DMatrix, DVector, MatrixMN, Dynamic, U3, Matrix1x3};
use core::num::FpCategory::Infinite;

// pub struct abosParameters {
//     zmin
// }


fn get_chebyshev_distance(x1: f64, x2: f64, y1: f64, y2: f64) -> f64 {
    if (x1 - x2).abs() > (y1 - y2).abs() {
        (x1 - x2).abs()
    } else {
        (y1 - y2).abs()
    }
}



//TODO make closure 
fn get_min_chebyshev_distance(xyz_points: &MatrixMN<f64, Dynamic, U3>) -> f64 {
    let minChebyshevDistance: f64 = INFINITY;
    let xyz_iter1 = xyz_points.row_iter();


    for (i, row) in xyz_iter1.enumerate(){
        let xyz_iter2 = xyz_points.row_iter();
        for (j,row2) in xyz_iter2.enumerate(){
            let test = row.clone_owned()-row2;

            println!("{}",test);
        }
    }

    return 10.0
}

pub struct ABOSInputs {
    points: Vec<na::Vector3<f64>>
}

pub struct ABOSGrid {
    filter: f64, //INPUT resolution parameter
    xyz_points: MatrixMN<f64, Dynamic, U3>, //INPUT all points XYZ
    x1: f64, //minx
    x2: f64, //maxx
    y1: f64, //miny
    y2: f64, //maxy
    z1:f64, //min z
    z2: f64, //max z
    // i1: i32, //xsize of grid
    // j1: i32, //ysize of grid
    // dx: f64, //size of grid on x
    // dy: f64, //size of grid on y
    // dmc: f64, //minimal chebyshev distance
    // P: DMatrix<f64>, //elements of the elevation matrix [i,j] size
    // DP: DMatrix<f64>, //auxiliary matrix same size as P
    // Z: DVector<f64>, //vector if z coordinates XYZ
    // DZ: DVector<f64>, //auxiliary vector same size as Z
    // NB: DMatrix<f64>, // Matrix of nearest points on grid. Containing indexes to nearest point in the XYZ array
    // K: DMatrix<f64>, // Grid distance of each grid to the point indexed in NB
    // kMax: f64,  //maximal element of matrix K
    // RS: f64, // Resolution of map
    //A->B means copy of A into B

}

pub fn compute_grid_dimensions(x1:f64,x2:f64,y1:f64,y2:f64,dmc:f64, filter: i32, points: &Vec<na::Vector3<f64>>) {
    //step 1: define which side is greater, by checking x extents vs y extents

    //step 2: grid size is defined as (length(LongerSide)/Dmc) where Dmc is the smallest gridded distance between points

    //step 3 find your grid size for your larger dimension set as il = i0*k largest possible while less than Filter

    //step 5 find your grid size for your smaller dimension such that it is as close to that of il as possible
}

impl ABOSGrid {
    pub fn new(points: Vec<Vec<f64>>, filter:f64) -> ABOSGrid {
        let mut vec = Vec::new(); //real function call would just pass the vector in here
        for x in points.iter() {
            for y in x.iter() {
                vec.push(*y);
            }
        }
        let point_count = points.len();
        let dm: DMatrix<_> = DMatrix::from_iterator(3, point_count, vec.into_iter());
        let fix_dm = dm.transpose();

        //step 1: make an array with all the points
        let mut xyz_points: MatrixMN<f64, Dynamic, U3> = na::convert(fix_dm);

        //step 2: get Point range information
        let (x1,x2,y1,y2,z1,z2) = get_ranges(&xyz_points);

        //step 3: get Chebyshev distance
        let dmc = get_min_chebyshev_distance(&xyz_points);
        //step 3: get the grid dimensions
        // let (i1,j1,dx,dy,minx,miny) = get_dimensions()

        ABOSGrid{
            filter, //INPUT resolution parameter
            xyz_points, //INPUT all points XYZ
            x1, //minx
            x2, //maxx
            y1, //miny
            y2, //maxy
            z1, //min z
            z2, //max z
            // i1, //xsize of grid
            // j1, //ysize of grid
            // dx, //size of grid on x
            // dy, //size of grid on y
            // dmc, //minimal chebyshev distance
            // P, //elements of the elevation matrix [i,j] size
            // DP, //auxiliary matrix same size as P
            // Z, //vector if z coordinates XYZ
            // DZ, //auxiliary vector same size as Z
            // NB, // Matrix of nearest points on grid. Containing indexes to nearest point in the XYZ array
            // K, // Grid distance of each grid to the point indexed in NB
            // kMax,  //maximal element of matrix K
            // RS, // Resolution of map
        }
    }
}


pub fn get_ranges(points: &MatrixMN<f64, Dynamic, U3> ) -> (f64,f64,f64,f64,f64,f64){
    let x1 = points.column(0).min();
    let x2 = points.column(0).max();
    let y1 = points.column(1).min();
    let y2 = points.column(1).max();
    let z1 = points.column(2).min();
    let z2 = points.column(2).max();

    return(x1,x2,y1,y2,z1,z2)
}

pub fn dynamicMatrixCheck() {
    // let (mut x1, mut x2, mut y1, mut y2, mut z1, mut z2) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //x & y min(1) and max(2)
    // let i1 = 4;
    // let j1 = 5;
    // let numPoints = 10;
    //
    // let mut P = DMatrix::from_element(i1, j1, 0.0);
    // let mut DP = DMatrix::from_element(i1, j1, 0.0);
    //
    // let Z = DVector::from_element(numPoints, 0.0);
    // let DZ = DVector::from_element(numPoints, 0.0);
    // let X = DVector::from_element(numPoints, 0.0);
    // let Y = DVector::from_element(numPoints, 0.0);
    //
    // let mut K = DMatrix::from_element(i1, j1, 0.0);
    // let mut kMax = K.max();
    // let mut NB = DMatrix::from_element(i1, j1, 0.0);
    // let mut filter = 1.0; // parameter
    // let mut RS = get_chebyshev_distance(x1, x2, y1, y2) / filter;
    //
    // x1 = X.min();
    // x2 = X.max();
    // y1 = Y.min();
    // y2 = Y.max();
    // z1 = Z.min();
    // z2 = Z.max();
    //
    // let dx = (x2 - x1) / ((i1 - 1) as f64);
    // let dy = (y2 - y1) / ((j1 - 1) as f64);
    //
    //
    // println!("RS {:?}", RS);
    // println!("x1 {:?}", x1);
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

