///Crate to generate a surface with the ABOS algorithm proposed in Art of Surface Interpolation byM gr. Miroslav Dressler, Ph.D.
/// http://m.dressler.sweb.cz/ABOS.htm
/// 

extern crate approx;
// For the macro relative_eq! in nalgebraImportCheck
extern crate nalgebra as na;

pub const INFINITY: f64 = 1.0f64 / 0.0f64;

use na::{DMatrix, DVector, MatrixMN, Dynamic, U3, Matrix1x3, U1, DimName, Dim, DefaultAllocator, Scalar};
use core::num::FpCategory::Infinite;
use na::base::allocator::Allocator;

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

pub struct ABOSInputs {
    points: Vec<na::Vector3<f64>>
}

pub struct ABOSGrid {
    degree: i8,
    filter: f64,
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
    i1: i32, //xsize of grid
    j1: i32, //ysize of grid
    dx: f64, //size of grid on x
    dy: f64, //size of grid on y
    // P: MatrixMN<f64,Dynamic,Dynamic>, //elements of the elevation matrix [i,j] size
    // DP: MatrixMN<f64,Dynamic,Dynamic>, //auxiliary matrix same size as P
    // Z: DVector<f64>, //vector if z coordinates XYZ
    // DZ: DVector<f64>, //auxiliary vector same size as Z
    // NBmin:i32,
    //     // NBmax:i32,
    // K: DMatrix<f64>, // Grid distance of each grid to the point indexed in NB
    // kMax: f64,  //maximal element of matrix K
    // RS: f64, // Resolution of map
    //A->B means copy of A into B

}

pub fn compute_grid_dimensions(x1: f64, x2: f64, y1: f64, y2: f64, dmc: f64, filter: f64)
-> (i32,i32,f64,f64){
    //step 1: Always assuming x side is greater

    //step 2: grid size is defined as i0=round  x21/ Dmc 

    let i0:i32 = f64::round((x2 - x1) / dmc) as i32;
    //step 3 find your grid size for your larger dimension set as i1 = i0*k largest possible while less than Filter
    let mut i1 :i32 = 0;
    let mut i = 0;
    loop {
        i+=1;
        let potential_val:i32 = i0 * i;
        if  filter > potential_val as f64 {
            i1 = potential_val;
        }else{
            break;
        }
    }
    //step 5 find your grid size for your smaller dimension such that it is as close to that of il as possible
    let j1 :i32 = f64::round((y2 - y1) / (x2 - x1) * (i1 as f64 - 1.0)) as i32;

    let dx = (x2-x1)/ i1 as f64;
    let dy = (y2-y1)/ j1 as f64;

    return (i1,j1,dx,dy)
}


// struct Container<T: Scalar, I1: Dim + DimName, J1: Dim + DimName>
//     where DefaultAllocator: Allocator<T, I1, J1>
// {
//     matrix: MatrixMN<T, I1, J1>
// }
//
// impl<T:Scalar,I1:Dim+DimName,J1:Dim+DimName> Container<T, I1, J1> {
//     pub fn new(defVal: &impl Scalar,sizeX:&impl Dim+DimName, sizeY:&impl Dim+DimName) -> Container<T,I1,J1> where
//     T: Scalar,
//     I1: Dim + DimName,
//     J1: Dim + DimName,
//
//     {
//
//     }
// }

impl ABOSGrid {
    pub fn new(points: Vec<Vec<f64>>, filter: f64, degree: i8) -> ABOSGrid {
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
        let (x1, x2, y1, y2, z1, z2) = get_ranges(&xyz_points);

        //step 3: get Chebyshev distance
        let dmc = get_min_chebyshev_distance(&xyz_points);
        //step 3: get the grid dimensions
        let (i1,j1,dx,dy) = compute_grid_dimensions(x1,x2,y1,y2,dmc,filter);

        //step 4: Create P and DP

        let P: MatrixMN<f64,Dynamic,Dynamic> = MatrixMN::from_element_generic(Dynamic::from_usize(i1 as usize), Dynamic::from_usize(j1 as usize), 0.0);

        // Pcontainer.matrix::

        println!("{:?}",P);
        // step 0: compute Q R and L
        ABOSGrid {
            degree,
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

